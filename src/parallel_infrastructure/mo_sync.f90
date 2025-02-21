! The module <i>mo_communication</i>
! provides functionality for boundary exchange and global sums.
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


!--------------------------------------------------
! timers definition
!needs:
!   USE mo_timer, ONLY: timer_start, timer_stop, timers_level, <timers_names>...
!



!----------------------------
MODULE mo_sync
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------


USE mo_kind,               ONLY: sp, wp, dp, i8
USE mo_exception,          ONLY: finish, message, message_text
USE mo_model_domain,       ONLY: t_patch
USE mo_decomposition_tools,ONLY: t_grid_domain_decomp_info, get_local_index, &
  &                              get_valid_local_index
USE mo_parallel_config,    ONLY: nproma
USE mo_run_config,         ONLY: msg_level
USE mo_math_constants,     ONLY: pi
USE mo_impl_constants,     ONLY: min_rlcell_int, min_rledge_int, min_rlvert_int, max_dom
USE mo_impl_constants_grf, ONLY: grf_bdywidth_c, grf_bdywidth_e
USE mo_io_units,           ONLY: find_next_free_unit, filename_max
USE mo_mpi,                ONLY: p_pe, p_bcast, p_sum, p_max, p_min, p_send, p_recv,               &
  &                              p_comm_work_test,  p_comm_work, p_n_work, my_process_is_mpi_test, &
  &                              get_my_mpi_all_id, process_mpi_all_test_id,                       &
  &                              my_process_is_mpi_parallel, p_work_pe0,p_pe_work,                 &
  &                              comm_lev, glob_comm, comm_proc0,   &
  &                              p_gather, p_gatherv, num_test_procs
USE mo_parallel_config, ONLY:p_test_run,   &
  & n_ghost_rows, l_log_checks, l_fast_sum
USE mo_communication,      ONLY: exchange_data, exchange_data_4de1,            &
                                 exchange_data_mult, t_comm_pattern,           &
                                 blk_no, idx_no, idx_1d, get_np_recv,          &
                                 get_np_send, get_pelist_recv, exchange_data_mult_mixprec

USE mo_timer,           ONLY: timer_start, timer_stop, activate_sync_timers, &
  & timer_global_sum, timer_omp_global_sum, timer_ordglb_sum!, timer_omp_ordglb_sum
USE mo_fortran_tools,   ONLY: t_ptr_3d, t_ptr_3d_sp, insert_dimension, assert_acc_host_only, set_acc_host_or_device

IMPLICIT NONE

PRIVATE


!modules interface-------------------------------------------
!subroutines
PUBLIC :: sync_patch_array, check_patch_array, sync_idx,              &
          global_sum_array, omp_global_sum_array,                     &
          global_sum_array2, global_sum_array3, global_sum,           &
          sync_patch_array_mult, global_min, global_max,              &
          sync_patch_array_4de1, decomposition_statistics,            &
          enable_sync_checks, disable_sync_checks,                    &
          cumulative_sync_patch_array, complete_cumulative_sync,      &
          sync_patch_array_mult_mp

!
!variables

INTEGER, PARAMETER, PUBLIC :: SYNC_C = 1
INTEGER, PARAMETER, PUBLIC :: SYNC_E = 2
INTEGER, PARAMETER, PUBLIC :: SYNC_V = 3
INTEGER, PARAMETER, PUBLIC :: SYNC_C1 = 4


INTERFACE sync_patch_array
  MODULE PROCEDURE sync_patch_array_r2
  MODULE PROCEDURE sync_patch_array_r3
  MODULE PROCEDURE sync_patch_array_s3
  MODULE PROCEDURE sync_patch_array_i2
  MODULE PROCEDURE sync_patch_array_i3
  MODULE PROCEDURE sync_patch_array_l2
  MODULE PROCEDURE sync_patch_array_l3
END INTERFACE

INTERFACE check_patch_array
  MODULE PROCEDURE check_patch_array_2
  MODULE PROCEDURE check_patch_array_3
  MODULE PROCEDURE check_patch_array_4
END INTERFACE

INTERFACE global_min
  MODULE PROCEDURE global_min_0d
  MODULE PROCEDURE global_min_1d
  MODULE PROCEDURE global_min_0di
END INTERFACE

INTERFACE global_max
  MODULE PROCEDURE global_max_0d
  MODULE PROCEDURE global_max_1d
  MODULE PROCEDURE global_max_0di
END INTERFACE

INTERFACE global_sum
  MODULE PROCEDURE global_sum_0d
  MODULE PROCEDURE global_sum_0di
  MODULE PROCEDURE global_sum_1d
  MODULE PROCEDURE global_sum_1di
END INTERFACE

INTERFACE global_sum_array
  MODULE PROCEDURE global_sum_array_0d
  MODULE PROCEDURE global_sum_array_0di
  MODULE PROCEDURE global_sum_array_1d
  MODULE PROCEDURE global_sum_array_2d
  MODULE PROCEDURE global_sum_array_3d
END INTERFACE

INTERFACE omp_global_sum_array
  MODULE PROCEDURE omp_global_sum_array_1d
  MODULE PROCEDURE omp_global_sum_array_2d
END INTERFACE


! Unit for logging sync errors
INTEGER, SAVE :: log_unit = -1

! Flag if sync checks are enabled when sync_patch_array et al is called
LOGICAL, SAVE :: do_sync_checks = .TRUE.

!-------------------------------------------------------------------------

!> Type definition for "cumulative syncs": These are boundary
!  exchanges for 3D cell-based fields, collecting as many fields as
!  possible before actually performing the sync.
TYPE t_cumulative_sync
  REAL(wp),      POINTER :: f3d(:,:,:)
  TYPE(t_patch), POINTER :: p_patch
END TYPE t_cumulative_sync

!> max. no. fields that can be handled by "sync_patch_array_mult"
INTEGER, PARAMETER :: MAX_CUMULATIVE_SYNC = 5
!> No. of pending cumulative sync; dim (typ,patch_id,fieldno)
INTEGER :: ncumul_sync(4,max_dom) = 0
!> List of cumulative sync fields:
TYPE(t_cumulative_sync) :: cumul_sync(4,max_dom,MAX_CUMULATIVE_SYNC)

  INTERFACE sync_patch_array_mult
    MODULE PROCEDURE sync_patch_array_mult_dp
    MODULE PROCEDURE sync_patch_array_mult_sp
  END INTERFACE sync_patch_array_mult

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sync'

CONTAINS


!-------------------------------------------------------------------------
!
!> enable_sync_checks: Enables sync checks when sync_patch_array et al is called

SUBROUTINE enable_sync_checks()
  do_sync_checks = .TRUE.
END SUBROUTINE enable_sync_checks
!-------------------------------------------------------------------------
!
!> disable_sync_checks: Disables sync checks when sync_patch_array et al is called

SUBROUTINE disable_sync_checks()
  do_sync_checks = .FALSE.
END SUBROUTINE disable_sync_checks

  FUNCTION comm_pat_of_type(p_patch, typ) RESULT(p_pat)
   INTEGER, INTENT(IN) :: typ
   TYPE(t_patch), TARGET, INTENT(IN) :: p_patch
   CLASS(t_comm_pattern), POINTER :: p_pat
   IF(typ == SYNC_C) THEN
     p_pat => p_patch%comm_pat_c
   ELSE IF(typ == SYNC_E) THEN
     p_pat => p_patch%comm_pat_e
   ELSE IF(typ == SYNC_V) THEN
     p_pat => p_patch%comm_pat_v
   ELSE IF(typ == SYNC_C1) THEN
     p_pat => p_patch%comm_pat_c1
   ELSE
     CALL finish('comm_pat_of_type','Illegal type parameter')
   ENDIF
  END FUNCTION comm_pat_of_type

!-------------------------------------------------------------------------
!> Does boundary exchange for a 3-D REAL array.
!
SUBROUTINE sync_patch_array_r3(typ, p_patch, arr, opt_varname)
   INTEGER,       INTENT(IN)    :: typ
   TYPE(t_patch), TARGET, INTENT(IN) :: p_patch
   REAL(wp),      INTENT(INOUT) :: arr(:,:,:)
   CHARACTER(len=*), TARGET, INTENT(IN), OPTIONAL :: opt_varname
   CLASS(t_comm_pattern), POINTER :: p_pat

   ! If this is a verification run, check consistency before doing boundary exchange
   IF (p_test_run .AND. do_sync_checks) &
     CALL check_patch_array_3(typ, p_patch, arr, opt_varname)

   ! Boundary exchange for work PEs
    IF(my_process_is_mpi_parallel()) THEN
      p_pat => comm_pat_of_type(p_patch, typ)
      CALL exchange_data(p_pat, arr)
    ENDIF
END SUBROUTINE sync_patch_array_r3


!-------------------------------------------------------------------------
!> Does boundary exchange for a 3-D single precision array.
!
SUBROUTINE sync_patch_array_s3(typ, p_patch, arr, opt_varname)
   INTEGER,       INTENT(IN)    :: typ
   TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch
   REAL(sp),      INTENT(INOUT) :: arr(:,:,:)
   CHARACTER(len=*), TARGET, INTENT(IN), OPTIONAL :: opt_varname
   CLASS(t_comm_pattern), POINTER :: p_pat
   CHARACTER(len=4), SAVE, TARGET :: default_varname = 'sync'
   CHARACTER(len=:), POINTER :: varname

   ! If this is a verification run, check consistency before doing boundary exchange
   IF (p_test_run .AND. do_sync_checks) THEN
     IF(PRESENT(opt_varname)) THEN
       varname => opt_varname
     ELSE
       varname => default_varname
     ENDIF
     CALL check_patch_array_sp(typ, p_patch, arr, varname)
   ENDIF

   ! Boundary exchange for work PEs
   IF(my_process_is_mpi_parallel()) THEN
      p_pat => comm_pat_of_type(p_patch, typ)
      CALL exchange_data(p_pat, arr)
   ENDIF
END SUBROUTINE sync_patch_array_s3


!-------------------------------------------------------------------------
!> Does boundary exchange for a 3-D INTEGER array.
!
!  @note This implementation does not perform a consistency check
!        (p_test_run)!
!
SUBROUTINE sync_patch_array_i3(typ, p_patch, arr)
   INTEGER,       INTENT(IN)    :: typ
   TYPE(t_patch), TARGET, INTENT(IN) :: p_patch
   INTEGER,       INTENT(INOUT) :: arr(:,:,:)
   CLASS(t_comm_pattern), POINTER :: p_pat

   ! Boundary exchange for work PEs
   IF(my_process_is_mpi_parallel()) THEN
      p_pat => comm_pat_of_type(p_patch, typ)
      CALL exchange_data(p_pat, arr)
   ENDIF
END SUBROUTINE sync_patch_array_i3

  SUBROUTINE sync_patch_array_l3(typ, p_patch, arr)
    INTEGER,       INTENT(IN)    :: typ
    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch
    LOGICAL,       INTENT(INOUT) :: arr(:,:,:)
    CLASS(t_comm_pattern), POINTER :: p_pat

    ! Boundary exchange for work PEs
    IF(my_process_is_mpi_parallel()) THEN
      p_pat => comm_pat_of_type(p_patch, typ)
      CALL exchange_data(p_pat, arr)
    ENDIF
  END SUBROUTINE sync_patch_array_l3


!-------------------------------------------------------------------------
!> Does boundary exchange for a 2-D REAL array.
!
SUBROUTINE sync_patch_array_r2(typ, p_patch, arr, opt_varname)
   INTEGER,       INTENT(IN)    :: typ
   TYPE(t_patch), INTENT(IN) :: p_patch
   REAL(wp), TARGET, INTENT(INOUT) :: arr(:,:)
   CHARACTER*(*), INTENT(IN), OPTIONAL :: opt_varname
   ! local variable
   REAL(wp), POINTER :: arr3(:,:,:)

   CALL insert_dimension(arr3, arr, 2)
   CALL sync_patch_array_r3(typ, p_patch, arr3, opt_varname)
END SUBROUTINE sync_patch_array_r2


!-------------------------------------------------------------------------
!> Does boundary exchange for a 2-D INTEGER array.
!
!  @note This implementation does not perform a consistency check
!        (p_test_run)!
!
SUBROUTINE sync_patch_array_i2(typ, p_patch, arr)
   INTEGER,       INTENT(IN)    :: typ
   TYPE(t_patch), INTENT(IN)    :: p_patch
   INTEGER, TARGET, INTENT(INOUT) :: arr(:,:)
   ! local variable
   INTEGER, POINTER :: arr3(:,:,:)

   CALL insert_dimension(arr3, arr, 2)
   CALL sync_patch_array_i3(typ, p_patch, arr3)
END SUBROUTINE sync_patch_array_i2

  SUBROUTINE sync_patch_array_l2(typ, p_patch, arr)
    INTEGER,       INTENT(IN)    :: typ
    TYPE(t_patch), INTENT(IN)    :: p_patch
    LOGICAL, TARGET, INTENT(INOUT) :: arr(:,:)
    ! local variable
    LOGICAL, POINTER :: arr3(:,:,:)

    CALL insert_dimension(arr3, arr, 2)
    CALL sync_patch_array_l3(typ, p_patch, arr3)
  END SUBROUTINE sync_patch_array_l2


!-------------------------------------------------------------------------
!! Does boundary exchange for up to 5 3D cell-based fields and/or a 4D field.
!! The 4D field can alternatively be passed as an array of 3D fields.
!!
SUBROUTINE sync_patch_array_mult_dp(typ, p_patch, nfields, f3din1, f3din2, f3din3, &
                                 f3din4, f3din5, f4din, f3din_arr, opt_varname)

   INTEGER, INTENT(IN)             :: typ
   TYPE(t_patch), INTENT(IN), TARGET :: p_patch
   INTEGER,     INTENT(IN)         :: nfields

   REAL(dp), OPTIONAL, INTENT(INOUT) ::  f3din1(:,:,:), f3din2(:,:,:), f3din3(:,:,:), &
      &                                  f3din4(:,:,:), f3din5(:,:,:), f4din(:,:,:,:)
   TYPE(t_ptr_3d), OPTIONAL, INTENT(INOUT) :: f3din_arr(:)

   CLASS(t_comm_pattern), POINTER :: p_pat
   CHARACTER(len=*), TARGET, INTENT(IN), OPTIONAL :: opt_varname
   INTEGER :: i
   INTEGER :: ndim2tot ! Sum of second dimensions over all input fields

!-----------------------------------------------------------------------

   p_pat => comm_pat_of_type(p_patch, typ)

   ! If this is a verification run, check consistency before doing boundary exchange
   IF (p_test_run .AND. do_sync_checks) THEN
     IF (PRESENT(f4din)) THEN
       DO i = 1, SIZE(f4din,4)
         CALL check_patch_array_3(typ, p_patch, f4din(:,:,:,i), opt_varname)
       ENDDO
     ENDIF

     IF (PRESENT(f3din_arr)) THEN
       DO i = 1, SIZE(f3din_arr)
         CALL check_patch_array_3(typ, p_patch, f3din_arr(i)%p, opt_varname)
       ENDDO
     ENDIF

     IF (PRESENT(f3din1)) CALL check_patch_array_3(typ, p_patch, f3din1, opt_varname)
     IF (PRESENT(f3din2)) CALL check_patch_array_3(typ, p_patch, f3din2, opt_varname)
     IF (PRESENT(f3din3)) CALL check_patch_array_3(typ, p_patch, f3din3, opt_varname)
     IF (PRESENT(f3din4)) CALL check_patch_array_3(typ, p_patch, f3din4, opt_varname)
     IF (PRESENT(f3din5)) CALL check_patch_array_3(typ, p_patch, f3din5, opt_varname)
   ENDIF

   ! Boundary exchange for work PEs
   IF(my_process_is_mpi_parallel()) THEN
     IF (PRESENT(f4din)) THEN
       ndim2tot = SIZE(f4din,4)*SIZE(f4din,2)
     ELSE
       ndim2tot = 0
     ENDIF
     IF (PRESENT(f3din_arr)) THEN
       DO i = 1, SIZE(f3din_arr)
         ndim2tot = ndim2tot + SIZE(f3din_arr(i)%p,2)
       ENDDO
     ENDIF
     IF (PRESENT(f3din1)) ndim2tot = ndim2tot+SIZE(f3din1,2)
     IF (PRESENT(f3din2)) ndim2tot = ndim2tot+SIZE(f3din2,2)
     IF (PRESENT(f3din3)) ndim2tot = ndim2tot+SIZE(f3din3,2)
     IF (PRESENT(f3din4)) ndim2tot = ndim2tot+SIZE(f3din4,2)
     IF (PRESENT(f3din5)) ndim2tot = ndim2tot+SIZE(f3din5,2)

     CALL exchange_data_mult(p_pat, nfields, ndim2tot, recv1=f3din1, recv2=f3din2, &
       &                     recv3=f3din3, recv4=f3din4, recv5=f3din5, recv4d=f4din , &
       &                     recv3d_arr=f3din_arr)
   ENDIF

END SUBROUTINE sync_patch_array_mult_dp

!-------------------------------------------------------------------------
!! Does boundary exchange for up to 5 3D cell-based fields and/or a 4D field.
!! The 4D field can alternatively be passed as an array of 3D fields.
!!
SUBROUTINE sync_patch_array_mult_sp(typ, p_patch, nfields, f3din1, f3din2, f3din3, &
                                 f3din4, f3din5, f4din, f3din_arr, opt_varname)

   INTEGER, INTENT(IN)             :: typ
   TYPE(t_patch), INTENT(IN), TARGET :: p_patch
   INTEGER,     INTENT(IN)         :: nfields

   REAL(sp), TARGET, INTENT(INOUT) :: f3din1(:,:,:)
   REAL(sp), TARGET, OPTIONAL, INTENT(INOUT) :: f3din2(:,:,:), f3din3(:,:,:), &
      &                           f3din4(:,:,:), f3din5(:,:,:), f4din(:,:,:,:)
   TYPE(t_ptr_3d_sp), OPTIONAL, INTENT(INOUT) :: f3din_arr(:)

   CLASS(t_comm_pattern), POINTER :: p_pat
   CHARACTER(len=*), TARGET, INTENT(IN), OPTIONAL :: opt_varname
   INTEGER :: i, nfields_
   INTEGER :: ndim2tot ! Sum of second dimensions over all input fields
   TYPE(t_ptr_3d_sp) :: fld(nfields)
   CHARACTER(len=*), PARAMETER :: routine &
        = modname//'::sync_patch_array_mult_sp'
!-----------------------------------------------------------------------

   p_pat => comm_pat_of_type(p_patch, typ)

   nfields_ = 0
   nfields_ = nfields_ + 1
   ndim2tot = SIZE(f3din1,2)
   fld(nfields_)%p => f3din1
   IF (PRESENT(f3din2)) THEN
     nfields_ = nfields_ + 1
     ndim2tot = ndim2tot + SIZE(f3din2,2)
     fld(nfields_)%p => f3din2
   END IF
   IF (PRESENT(f3din3)) THEN
     nfields_ = nfields_ + 1
     ndim2tot = ndim2tot + SIZE(f3din3,2)
     fld(nfields_)%p => f3din3
   END IF
   IF (PRESENT(f3din4)) THEN
     nfields_ = nfields_ + 1
     ndim2tot = ndim2tot + SIZE(f3din4,2)
     fld(nfields_)%p => f3din4
   END IF
   IF (PRESENT(f3din5)) THEN
     nfields_ = nfields_ + 1
     ndim2tot = ndim2tot + SIZE(f3din5,2)
     fld(nfields_)%p => f3din5
   END IF
   IF (PRESENT(f4din)) THEN
     DO i = 1, SIZE(f4din,4)
       nfields_ = nfields_ + 1
       fld(nfields_)%p => f4din(:,:,:,i)
     ENDDO
   ENDIF
   IF (PRESENT(f3din_arr)) THEN
     DO i = 1, SIZE(f3din_arr)
       ndim2tot = ndim2tot + SIZE(f3din_arr(i)%p, 2)
       fld(nfields_+i)%p => f3din_arr(i)%p
     ENDDO
   ENDIF
   IF (nfields_ /= nfields) THEN
     CALL finish(routine, 'internal error')
   END IF
   ! If this is a verification run, check consistency before doing boundary exchange
   IF (p_test_run .AND. do_sync_checks) THEN
     DO i = 1, nfields
       CALL check_patch_array_sp(typ, p_patch, fld(i)%p, opt_varname)
     ENDDO
   ENDIF

   ! Boundary exchange for work PEs
   IF(my_process_is_mpi_parallel()) THEN
     CALL p_pat%exchange_data_mult_mixprec(0, 0, nfields, ndim2tot, &
          recv_sp=fld)
   ENDIF

 END SUBROUTINE sync_patch_array_mult_sp


!-------------------------------------------------------------------------
!! Does boundary exchange for up to 5 3D cell-based fields and/or a 4D field,
!! which can either be single precision or double precision
!!
SUBROUTINE sync_patch_array_mult_mp(typ, p_patch, nfields, nfields_sp, f3din1, f3din2, f3din3, &
  f3din4, f3din5, f3din1_sp, f3din2_sp, f3din3_sp, f3din4_sp, f3din5_sp, f4din, f4din_sp, opt_varname)

   INTEGER, INTENT(IN)               :: typ
   TYPE(t_patch), INTENT(IN), TARGET :: p_patch
   INTEGER,     INTENT(IN)           :: nfields, nfields_sp

   REAL(dp), OPTIONAL, INTENT(INOUT) ::  f3din1(:,:,:), f3din2(:,:,:), f3din3(:,:,:), &
      &                                  f3din4(:,:,:), f3din5(:,:,:), f4din(:,:,:,:)
   REAL(sp), OPTIONAL, INTENT(INOUT) ::  f3din1_sp(:,:,:), f3din2_sp(:,:,:), f3din3_sp(:,:,:), &
      &                                  f3din4_sp(:,:,:), f3din5_sp(:,:,:), f4din_sp(:,:,:,:)

   CLASS(t_comm_pattern), POINTER :: p_pat
   CHARACTER(len=*), TARGET, INTENT(IN), OPTIONAL :: opt_varname

   INTEGER :: i
   INTEGER :: ndim2tot, ndim2tot_sp ! Sum of second dimensions over all input fields

   CHARACTER(len=4), SAVE, TARGET :: default_name = 'sync'

!-----------------------------------------------------------------------

    IF(typ == SYNC_C) THEN
      p_pat => p_patch%comm_pat_c
    ELSE IF(typ == SYNC_E) THEN
      p_pat => p_patch%comm_pat_e
    ELSE IF(typ == SYNC_V) THEN
      p_pat => p_patch%comm_pat_v
    ELSE IF(typ == SYNC_C1) THEN
      p_pat => p_patch%comm_pat_c1
    ENDIF

   ! If this is a verification run, check consistency before doing boundary exchange
   IF (p_test_run .AND. do_sync_checks) THEN
     IF (PRESENT(f4din)) THEN
       DO i = 1, SIZE(f4din,4)
         CALL check_patch_array_3(typ, p_patch, f4din(:,:,:,i), opt_varname)
       ENDDO
     ENDIF
     IF (PRESENT(f4din_sp)) THEN
       DO i = 1, SIZE(f4din_sp,4)
         CALL check_patch_array_sp(typ, p_patch, f4din_sp(:,:,:,i), opt_varname)
       ENDDO
     ENDIF
     IF (PRESENT(f3din1)) CALL check_patch_array_3(typ, p_patch, f3din1, opt_varname)
     IF (PRESENT(f3din2)) CALL check_patch_array_3(typ, p_patch, f3din2, opt_varname)
     IF (PRESENT(f3din3)) CALL check_patch_array_3(typ, p_patch, f3din3, opt_varname)
     IF (PRESENT(f3din4)) CALL check_patch_array_3(typ, p_patch, f3din4, opt_varname)
     IF (PRESENT(f3din5)) CALL check_patch_array_3(typ, p_patch, f3din5, opt_varname)
     IF (PRESENT(f3din1_sp)) CALL check_patch_array_sp(typ, p_patch, f3din1_sp, opt_varname)
     IF (PRESENT(f3din2_sp)) CALL check_patch_array_sp(typ, p_patch, f3din2_sp, opt_varname)
     IF (PRESENT(f3din3_sp)) CALL check_patch_array_sp(typ, p_patch, f3din3_sp, opt_varname)
     IF (PRESENT(f3din4_sp)) CALL check_patch_array_sp(typ, p_patch, f3din4_sp, opt_varname)
     IF (PRESENT(f3din5_sp)) CALL check_patch_array_sp(typ, p_patch, f3din5_sp, opt_varname)
   ENDIF

   ! Boundary exchange for work PEs
   IF(my_process_is_mpi_parallel()) THEN
     IF (PRESENT(f4din)) THEN
       ndim2tot = SIZE(f4din,4)*SIZE(f4din,2)
     ELSE
       ndim2tot = 0
     ENDIF
     IF (PRESENT(f3din1)) ndim2tot = ndim2tot+SIZE(f3din1,2)
     IF (PRESENT(f3din2)) ndim2tot = ndim2tot+SIZE(f3din2,2)
     IF (PRESENT(f3din3)) ndim2tot = ndim2tot+SIZE(f3din3,2)
     IF (PRESENT(f3din4)) ndim2tot = ndim2tot+SIZE(f3din4,2)
     IF (PRESENT(f3din5)) ndim2tot = ndim2tot+SIZE(f3din5,2)

     IF (PRESENT(f4din_sp)) THEN
       ndim2tot_sp = SIZE(f4din_sp,4)*SIZE(f4din_sp,2)
     ELSE
       ndim2tot_sp = 0
     ENDIF
     IF (PRESENT(f3din1_sp)) ndim2tot_sp = ndim2tot_sp+SIZE(f3din1_sp,2)
     IF (PRESENT(f3din2_sp)) ndim2tot_sp = ndim2tot_sp+SIZE(f3din2_sp,2)
     IF (PRESENT(f3din3_sp)) ndim2tot_sp = ndim2tot_sp+SIZE(f3din3_sp,2)
     IF (PRESENT(f3din4_sp)) ndim2tot_sp = ndim2tot_sp+SIZE(f3din4_sp,2)
     IF (PRESENT(f3din5_sp)) ndim2tot_sp = ndim2tot_sp+SIZE(f3din5_sp,2)

     CALL exchange_data_mult_mixprec(p_pat, nfields, ndim2tot, nfields_sp, ndim2tot_sp,                    &
       recv1_dp=f3din1,    recv2_dp=f3din2,    recv3_dp=f3din3,    recv4_dp=f3din4,    recv5_dp=f3din5,    &
       recv1_sp=f3din1_sp, recv2_sp=f3din2_sp, recv3_sp=f3din3_sp, recv4_sp=f3din4_sp, recv5_sp=f3din5_sp, &
       recv4d_dp=f4din,    recv4d_sp=f4din_sp                                                              )
   ENDIF

END SUBROUTINE sync_patch_array_mult_mp


!! Does boundary exchange for a 4D field for which the extra dimension
!! is on the third index.
!!
SUBROUTINE sync_patch_array_4de1(typ, p_patch, nfields, f4din, opt_varname)

   INTEGER, INTENT(IN)             :: typ
   TYPE(t_patch), INTENT(IN), TARGET :: p_patch
   INTEGER,     INTENT(IN)         :: nfields

   REAL(wp), INTENT(INOUT) :: f4din(:,:,:,:)

   CLASS(t_comm_pattern), POINTER :: p_pat
   CHARACTER(len=*), TARGET, INTENT(IN), OPTIONAL :: opt_varname
   INTEGER :: i, ndim2tot

!-----------------------------------------------------------------------

    p_pat => comm_pat_of_type(p_patch, typ)

   ! If this is a verification run, check consistency before doing boundary exchange
   IF (p_test_run .AND. do_sync_checks) THEN
! The silly PGI OpenACC compiler does not know f4din(i,:,:,:) is present on the device
     DO i = 1, nfields
       CALL check_patch_array_3(typ, p_patch, f4din(i,:,:,:), opt_varname)
     ENDDO
   ENDIF

   ! Boundary exchange for work PEs
     IF(my_process_is_mpi_parallel()) THEN
     IF (nfields/=UBOUND(f4din,1)) &
       CALL finish('sync_patch_array_4de1','inconsistent arguments')
     ndim2tot = nfields*SIZE(f4din,3)
     CALL exchange_data_4de1(p_pat, nfields, ndim2tot, recv=f4din)
   ENDIF

END SUBROUTINE sync_patch_array_4de1



!-------------------------------------------------------------------------
!
!

!! Wrapper routine for checking single precision arrays
SUBROUTINE check_patch_array_sp(typ, p_patch, arr, opt_varname)

   INTEGER, INTENT(IN)     :: typ
   TYPE(t_patch), INTENT(IN), TARGET :: p_patch
   CHARACTER(*), INTENT(IN), OPTIONAL :: opt_varname

   REAL(sp), INTENT(IN) :: arr(:,:,:)
   REAL(wp) :: arr_wp(SIZE(arr,1),SIZE(arr,2),SIZE(arr,3))

   !$ACC DATA CREATE(arr_wp) IF(i_am_accel_node)
   !$ACC KERNELS ASYNC(1) IF(i_am_accel_node)
   arr_wp(:,:,:) = REAL(arr(:,:,:),wp)
   !$ACC END KERNELS
   CALL check_patch_array_3(typ, p_patch, arr_wp, opt_varname)
   !$ACC WAIT(1)
   !$ACC END DATA

END SUBROUTINE check_patch_array_sp


!! In a verification run, this routine checks the consistency of an array,
!! i.e. if the parts on the worker PEs are identical with the data on
!! the verification PE.
!! For a non-verification run it just does nothing.
!!
SUBROUTINE check_patch_array_3(typ, p_patch, arr, opt_varname)

!

   INTEGER, INTENT(IN)     :: typ
   TYPE(t_patch), INTENT(IN), TARGET :: p_patch
   CHARACTER(*), INTENT(IN), OPTIONAL :: opt_varname

   REAL(wp), INTENT(IN) :: arr(:,:,:)

   REAL(wp), ALLOCATABLE:: arr_g(:,:,:)
   INTEGER :: j, jb, jl, jb_g, jl_g, n, ndim2, ndim3, nblks_g, flag, jk
   INTEGER :: ityp, ndim, ndim_g, jk_min_err
   INTEGER :: nerr(0:n_ghost_rows), shape_recv(3)
   INTEGER, POINTER :: p_glb_index(:), p_decomp_domain(:,:)
   CLASS(t_comm_pattern), POINTER :: p_pat_work2test
   LOGICAL :: l_my_process_is_mpi_test

   CHARACTER(len=256) :: varname, cfmt
   INTEGER :: varname_tlen

   CHARACTER(filename_max) :: log_file
   REAL(wp) :: absmax, relmax
   LOGICAL :: sync_error

   ityp   = -1
   ndim   = -1
   ndim_g = -1
   sync_error = .FALSE.

   NULLIFY(p_glb_index, p_decomp_domain)
!-----------------------------------------------------------------------

   IF(.NOT. p_test_run) RETURN ! This routine is only effective in a verification run

   IF(PRESENT(opt_varname)) THEN
     varname = opt_varname
     varname_tlen = LEN(opt_varname)
   ELSE
     varname = ' no VARNAME supplied'
     varname_tlen = 20
   ENDIF

   ! Check dimensions of arr, determine if this is an cell/edge/vert array

   IF(UBOUND(arr,1) /= nproma) THEN
      CALL finish('sync_patch_array','first dimension /= nproma')
   ENDIF

   ndim2 = UBOUND(arr,2)
   ndim3 = UBOUND(arr,3)

   !$ACC DATA PRESENT(arr) IF(i_am_accel_node)
   !$ACC UPDATE HOST(arr) ASYNC(1) IF(i_am_accel_node)

   IF(typ == SYNC_C .OR. typ == SYNC_C1) THEN
      ndim   = p_patch%n_patch_cells
      ndim_g = p_patch%n_patch_cells_g
      p_glb_index => p_patch%cells%decomp_info%glb_index
      p_decomp_domain => p_patch%cells%decomp_info%decomp_domain
      ityp = typ
      p_pat_work2test => p_patch%comm_pat_work2test(1)%p
   ELSE IF(typ == SYNC_E) THEN
      ndim   = p_patch%n_patch_edges
      ndim_g = p_patch%n_patch_edges_g
      p_glb_index => p_patch%edges%decomp_info%glb_index
      p_decomp_domain => p_patch%edges%decomp_info%decomp_domain
      ityp = typ
      p_pat_work2test => p_patch%comm_pat_work2test(3)%p
   ELSE IF(typ == SYNC_V) THEN
      ndim   = p_patch%n_patch_verts
      ndim_g = p_patch%n_patch_verts_g
      p_glb_index => p_patch%verts%decomp_info%glb_index
      p_decomp_domain => p_patch%verts%decomp_info%decomp_domain
      ityp = typ
      p_pat_work2test => p_patch%comm_pat_work2test(2)%p
   ELSE IF(typ == 0) THEN
      ! typ == 0 may be set for quick checks without knowing the type of the array.
      ! It may only be used if the array is correctly dimensioned.

      IF(ndim3 == p_patch%nblks_c) THEN
         ndim   = p_patch%n_patch_cells
         ndim_g = p_patch%n_patch_cells_g
         p_glb_index => p_patch%cells%decomp_info%glb_index
         p_decomp_domain => p_patch%cells%decomp_info%decomp_domain
         ityp = SYNC_C
         p_pat_work2test => p_patch%comm_pat_work2test(1)%p
      ELSE IF(ndim3 == p_patch%nblks_e) THEN
         ndim   = p_patch%n_patch_edges
         ndim_g = p_patch%n_patch_edges_g
         p_glb_index => p_patch%edges%decomp_info%glb_index
         p_decomp_domain => p_patch%edges%decomp_info%decomp_domain
         ityp = SYNC_E
         p_pat_work2test => p_patch%comm_pat_work2test(3)%p
      ELSE IF(ndim3 == p_patch%nblks_v) THEN
         ndim   = p_patch%n_patch_verts
         ndim_g = p_patch%n_patch_verts_g
         p_glb_index => p_patch%verts%decomp_info%glb_index
         p_decomp_domain => p_patch%verts%decomp_info%decomp_domain
         ityp = SYNC_V
         p_pat_work2test => p_patch%comm_pat_work2test(2)%p
      ELSE
         CALL finish('check_patch_array','typ==0 but unknown blocksize of array')
      ENDIF
   ELSE
      CALL finish('sync_patch_array','Illegal type parameter')
   ENDIF

   ! Actually do the check.
   ! The test PE broadcasts its full array, the other check if their section matches.

   nblks_g = (ndim_g-1)/nproma+1

   l_my_process_is_mpi_test = my_process_is_mpi_test()
   !$ACC WAIT(1) !GV: UPDATE HOST(arr) finished
   IF (num_test_procs > 1) THEN
     shape_recv = SHAPE(arr)
     ALLOCATE(arr_g(shape_recv(1),shape_recv(2),shape_recv(3)))
     CALL exchange_data(p_pat_work2test, arr_g, arr)
     IF(l_my_process_is_mpi_test) THEN
       jk_min_err = HUGE(jk_min_err)
!$OMP PARALLEL PRIVATE(jb,jk,jl) REDUCTION(.or.: sync_error) &
!$OMP REDUCTION(MIN: jk_min_err)
!$OMP DO
       DO jb = 1, ndim3
         DO jk = 1, ndim2
           DO jl = 1, nproma
             IF (p_decomp_domain(jl,jb) == 0) THEN
               sync_error = sync_error &
                    .OR. arr(jl, jk, jb) /= arr_g(jl, jk, jb)
               jk_min_err = MIN(jk_min_err, MERGE(jk, jk_min_err, &
                    arr(jl, jk, jb) /= arr_g(jl, jk, jb)))
             END IF
           END DO
         END DO
       END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
     END IF
   ELSE IF(l_my_process_is_mpi_test) THEN

     ! the test PE may also have reordered global indices. create a
     ! temporary array in the correct order:
     ALLOCATE(arr_g(nproma,ndim2,nblks_g))
     DO j=1,ndim
         jb = blk_no(j) ! Block index in distributed patch
         jl = idx_no(j) ! Line  index in distributed patch

         jb_g = blk_no(p_glb_index(j)) ! Block index in global patch
         jl_g = idx_no(p_glb_index(j)) ! Line  index in global patch
         arr_g(jl_g,1:ndim2,jb_g) = arr(jl,1:ndim2,jb)
     END DO

     IF(comm_lev==0) THEN
       CALL p_bcast(arr_g(:,:,1:nblks_g), process_mpi_all_test_id, comm=p_comm_work_test)
     ELSE
       CALL p_send(arr_g(:,:,1:nblks_g),comm_proc0(comm_lev)+p_work_pe0,1)
     ENDIF
     DEALLOCATE(arr_g)

   ELSE

      ALLOCATE(arr_g(nproma,ndim2,nblks_g))
      IF(comm_lev==0) THEN
         CALL p_bcast(arr_g(:,:,1:nblks_g), process_mpi_all_test_id, comm=p_comm_work_test)
      ELSE
         IF(p_pe_work==comm_proc0(comm_lev)) &
           & CALL p_recv(arr_g(:,:,1:nblks_g), process_mpi_all_test_id, 1)
         CALL p_bcast(arr_g(:,:,1:nblks_g),0,comm=glob_comm(comm_lev))
      ENDIF

      ! Count errors in the inner domain and the different ghost rows

      nerr(:) = 0
      absmax = 0.0_wp
      relmax = 0.0_wp

      DO j = 1, ndim

         jb = blk_no(j) ! Block index in distributed patch
         jl = idx_no(j) ! Line  index in distributed patch

         jb_g = blk_no(p_glb_index(j)) ! Block index in global patch
         jl_g = idx_no(p_glb_index(j)) ! Line  index in global patch

         flag = p_decomp_domain(jl,jb)

         ! Safety measures only:
         flag = MAX(flag,0)
         flag = MIN(flag,UBOUND(nerr,1))

         DO n=1,ndim2
            IF(arr(jl,n,jb) /= arr_g(jl_g,n,jb_g)) THEN
               nerr(flag) = nerr(flag)+1
               IF(flag==0) THEN
                  ! Real sync error detected
                  sync_error = .TRUE.
                  absmax = MAX(absmax,ABS(arr(jl,n,jb) - arr_g(jl_g,n,jb_g)))
                  IF (l_log_checks) THEN





                     WRITE(log_unit,'(2a,5i7,3e18.10)') varname, 'sync error location:',&
                       jb,jl,jb_g,jl_g,n,arr(jl,n,jb),arr_g(jl_g,n,jb_g),    &
                       ABS(arr(jl,n,jb)-arr_g(jl_g,n,jb_g))

                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      IF(l_log_checks) THEN

         IF(log_unit<0) THEN
            WRITE(log_file,'(''log'',i4.4,''.txt'')') p_pe
            log_unit = find_next_free_unit(10,99)
            OPEN(log_unit,FILE=log_file)
         ENDIF

         n = n_ghost_rows
         WRITE(cfmt,'(a,i3,a)') '(',n+1,'i8,'' '',2a)'

         IF(ALL(arr == 0.0_wp)) THEN
            WRITE(log_unit,cfmt) nerr(0:n),varname(1:varname_tlen), ': ALL 0 !!!'
         ELSE
            WRITE(log_unit,cfmt) nerr(0:n),varname(1:varname_tlen)
         ENDIF



         IF(absmax > 0.0_wp) WRITE(log_unit,*) 'Max abs inner err:',absmax

      ENDIF

      ! Terminate the programm if the array is out of sync
      DEALLOCATE(arr_g)

   ENDIF
   IF (sync_error) THEN
     IF (num_test_procs > 1) &
          WRITE(0, '(2a,i0)') varname(1:varname_tlen), &
          ' sync error in level jk = ', jk_min_err
     IF(l_log_checks) THEN
       CLOSE (log_unit)
     ENDIF
     CALL finish('sync_patch_array','Out of sync detected!')
   ENDIF

   !$ACC END DATA

END SUBROUTINE check_patch_array_3
!-------------------------------------------------------------------------
!
!

!! 2-D Interface to check_patch_array.
!!
SUBROUTINE check_patch_array_2(typ, p_patch, arr, opt_varname)

!

   INTEGER, INTENT(IN)     :: typ
   TYPE(t_patch), INTENT(IN) :: p_patch
   REAL(wp), TARGET, INTENT(IN)    :: arr(:,:)
   CHARACTER(*), INTENT(IN), OPTIONAL :: opt_varname

   REAL(wp), POINTER :: arr3(:,:,:)
!-----------------------------------------------------------------------

   IF(.NOT. p_test_run) RETURN ! This routine is only effective in a verification run
   CALL insert_dimension(arr3, arr, 2)

   CALL check_patch_array_3(typ, p_patch, arr3, opt_varname)

END SUBROUTINE check_patch_array_2
!-------------------------------------------------------------------------
!
!

!! 4-D Interface to check_patch_array.
!!
SUBROUTINE check_patch_array_4(typ, p_patch, arr, opt_varname)

!

   INTEGER, INTENT(IN)     :: typ
   TYPE(t_patch), INTENT(IN), TARGET :: p_patch
   REAL(wp), INTENT(INOUT) :: arr(:,:,:,:)
   CHARACTER(len=*), INTENT(IN), OPTIONAL :: opt_varname

   CHARACTER(len=256) :: new_var
   INTEGER :: jt
!-----------------------------------------------------------------------

   IF(.NOT. p_test_run) RETURN ! This routine is only effective in a verification run

   DO jt=1,UBOUND(arr,4)
      IF(PRESENT(opt_varname)) THEN
         WRITE(new_var,'(a,''['',i2,'']'')') opt_varname, jt
         CALL check_patch_array_3(typ, p_patch, arr(:,:,:,jt),TRIM(new_var))
      ELSE
         CALL check_patch_array_3(typ, p_patch, arr(:,:,:,jt))
      ENDIF
   ENDDO

END SUBROUTINE check_patch_array_4
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!! Syncs an idx/blk pair of arrays

SUBROUTINE sync_idx(type_arr, type_idx, p_patch, idx, blk, opt_remap, opt_varname )

  INTEGER, INTENT(IN) :: type_arr, type_idx
  TYPE(t_patch), TARGET, INTENT(INOUT) :: p_patch
  INTEGER, INTENT(INOUT) :: idx(:,:), blk(:,:)
  LOGICAL, INTENT(IN), OPTIONAL :: opt_remap
  CHARACTER(len=*), TARGET, INTENT(IN), OPTIONAL :: opt_varname

  INTEGER :: nblks, n_idx, n_idx_g, jb, jl, i_l, i_g
  LOGICAL :: remap
  REAL(wp), ALLOCATABLE :: z_idx(:,:)
  TYPE(t_grid_domain_decomp_info), POINTER :: decomp_info


  ! opt_remap: Flag if index values pointing outside local domain
  ! should be remapped to values within local domain

  IF(PRESENT(opt_remap)) THEN
    remap = opt_remap
  ELSE
    remap = .TRUE.
  ENDIF

  IF(type_arr == SYNC_C) THEN
    nblks = p_patch%nblks_c
  ELSEIF(type_arr == SYNC_E) THEN
    nblks = p_patch%nblks_e
  ELSEIF(type_arr == SYNC_V) THEN
    nblks = p_patch%nblks_v
  ELSE
    CALL finish('sync_idx','Unsupported type_arr')
  ENDIF

  IF(type_idx == SYNC_C) THEN
    decomp_info => p_patch%cells%decomp_info
    n_idx = p_patch%n_patch_cells
    n_idx_g = p_patch%n_patch_cells_g
  ELSEIF(type_idx == SYNC_E) THEN
    decomp_info => p_patch%edges%decomp_info
    n_idx = p_patch%n_patch_edges
    n_idx_g = p_patch%n_patch_edges_g
  ELSEIF(type_idx == SYNC_V) THEN
    decomp_info => p_patch%verts%decomp_info
    n_idx = p_patch%n_patch_verts
    n_idx_g = p_patch%n_patch_verts_g
  ELSE
    CALL finish('sync_idx','Unsupported type_idx')
  ENDIF


  if(ubound(idx,1) /= nproma) CALL finish('sync_idx','ubound(idx,1) /= nproma')
  if(ubound(idx,2) /= nblks)  CALL finish('sync_idx','ubound(idx,2) /= nblks')

  ALLOCATE(z_idx(nproma,nblks))
  z_idx = 0._wp

  ! Set z_idx with the global 1D-index of all points

  !$ACC DATA COPYIN(z_idx) IF(i_am_accel_node)

  !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  DO jb = 1, nblks
    DO jl = 1, nproma

      i_l = idx_1d(idx(jl,jb),blk(jl,jb))

      IF(i_l <= 0 .or. i_l > n_idx) THEN
        z_idx(jl,jb) = 0._wp
      ELSE
        z_idx(jl,jb) = REAL(decomp_info%glb_index(i_l),wp)
      ENDIF

    END DO
  END DO
  !$ACC END PARALLEL
  !$ACC WAIT(1)

  ! Sync z_idx
  CALL sync_patch_array(type_arr, p_patch, z_idx, opt_varname)

  ! Set all points with local index corresponding to z_idx
  !$ACC PARALLEL ASYNC(1) IF(i_am_accel_node)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  DO jb = 1, nblks
    DO jl = 1, nproma

      i_g = INT(z_idx(jl,jb))

      IF(i_g <= 0 .or. i_g > n_idx_g) THEN
        idx(jl,jb) = idx_no(0)
        blk(jl,jb) = blk_no(0)
      ELSE
        IF (remap) THEN
          i_l = MAX(get_valid_local_index(decomp_info%glb2loc_index, i_g), 1)
        ELSE
          i_l = get_local_index(decomp_info%glb2loc_index, i_g)
          ! Set it to negative global index like in
          ! mo_setup_subdivision/get_local_idx_blk in case the local
          ! index is invalid
          ! MoHa: ...does this make sense?
          IF (i_l <= 0) i_l = - i_g
        END IF
        idx(jl,jb) = idx_no(i_l)
        blk(jl,jb) = blk_no(i_l)
      ENDIF

    END DO


  END DO
  !$ACC END PARALLEL
  !$ACC WAIT(1)

  !$ACC END DATA

END SUBROUTINE sync_idx

! Routines for global summation of scalar quantities or arrays of scalar quantities
! Unlike the subsequent global_sum_array routines, no summation is made over the array elements
!
FUNCTION global_sum_0d (z_in, opt_iroot) RESULT (global_sum)

  REAL(wp),          INTENT(in) :: z_in
  INTEGER, OPTIONAL,INTENT(IN)  :: opt_iroot
  REAL(wp)                      :: global_sum

  INTEGER :: p_comm_glob
!-----------------------------------------------------------------------

  IF(comm_lev==0) THEN
    p_comm_glob = p_comm_work
  ELSE
    p_comm_glob = glob_comm(comm_lev)
  ENDIF

  global_sum = p_sum(z_in, comm=p_comm_glob, root=opt_iroot)


END FUNCTION global_sum_0d

! integer variant of global_sum_0d
!
FUNCTION global_sum_0di (z_in, opt_iroot) RESULT (global_sum)

  INTEGER,          INTENT(in) :: z_in
  INTEGER, OPTIONAL,INTENT(IN) :: opt_iroot
  INTEGER                      :: global_sum
 
  INTEGER :: p_comm_glob, i_in(1), i_out(1)
!-----------------------------------------------------------------------

  IF(comm_lev==0) THEN
    p_comm_glob = p_comm_work
  ELSE
    p_comm_glob = glob_comm(comm_lev)
  ENDIF

  i_in(1)    = z_in
  i_out      = p_sum(i_in, comm=p_comm_glob, root=opt_iroot)
  global_sum = i_out(1)

END FUNCTION global_sum_0di

! Variant for 1D arrays
!
FUNCTION global_sum_1d (zfield, opt_iroot) RESULT (global_sum)

  REAL(wp),          INTENT(in) :: zfield(:)
  INTEGER, OPTIONAL,INTENT(IN)  :: opt_iroot
  REAL(wp)                      :: global_sum(SIZE(zfield))

  INTEGER :: p_comm_glob
!-----------------------------------------------------------------------

  IF(comm_lev==0) THEN
    p_comm_glob = p_comm_work
  ELSE
    p_comm_glob = glob_comm(comm_lev)
  ENDIF

  global_sum = p_sum(zfield, comm=p_comm_glob, root=opt_iroot)


END FUNCTION global_sum_1d


! integer variant of global_sum_1d
!
FUNCTION global_sum_1di (zfield, opt_iroot) RESULT (global_sum)

  INTEGER,          INTENT(in) :: zfield(:)
  INTEGER, OPTIONAL,INTENT(IN) :: opt_iroot
  INTEGER                      :: global_sum(SIZE(zfield))
 
  INTEGER :: p_comm_glob
!-----------------------------------------------------------------------

  IF(comm_lev==0) THEN
    p_comm_glob = p_comm_work
  ELSE
    p_comm_glob = glob_comm(comm_lev)
  ENDIF

  global_sum = p_sum(zfield, comm=p_comm_glob, root=opt_iroot)


END FUNCTION global_sum_1di

!-------------------------------------------------------------------------
!! Calculates the global sum of an integer scalar.
!! This routine should be called outside an OMP parallel Region!
!!
!! @param[in] opt_iroot (Optional:) root PE, otherwise we perform an
!!            ALL-TO-ALL operation.
!!
FUNCTION global_sum_array_0di (zfield, opt_iroot) RESULT (global_sum)

  INTEGER,           INTENT(in) :: zfield
  INTEGER, OPTIONAL, INTENT(IN) :: opt_iroot
  INTEGER                       :: global_sum
  ! local variables
  REAL(wp)                      :: z_aux, z_auxs
  INTEGER                       :: p_comm_glob
  IF (activate_sync_timers) CALL timer_start(timer_global_sum)

  IF(comm_lev==0) THEN
    p_comm_glob = p_comm_work
  ELSE
    p_comm_glob = glob_comm(comm_lev)
  ENDIF

  z_aux  =  REAL(zfield,wp)
  z_auxs = 0._wp
  z_auxs = p_sum(z_aux, comm=p_comm_glob, root=opt_iroot)

  global_sum = NINT(z_auxs)

  IF (activate_sync_timers) CALL timer_stop(timer_global_sum)
END FUNCTION global_sum_array_0di

!-------------------------------------------------------------------------
!! Calculates the global sum of zfield and checks for consistency
!! when doing a verification run.
!! This routine should be called outside an OMP parallel Region!
!!
!! @param[in] opt_iroot (Optional:) root PE, otherwise we perform an
!!            ALL-TO-ALL operation (for "l_fast_sum == .TRUE.").
!!
FUNCTION global_sum_array_0d (zfield, opt_iroot) RESULT (global_sum)

  REAL(wp), INTENT(IN)           :: zfield
  INTEGER,  INTENT(IN), OPTIONAL :: opt_iroot
  REAL(wp)                       :: global_sum
  ! local variables
  REAL(wp)                      :: sum_on_testpe(1), z_aux(1)
  INTEGER :: p_comm_glob

  IF(comm_lev==0) THEN
    p_comm_glob = p_comm_work
  ELSE
    p_comm_glob = glob_comm(comm_lev)
  ENDIF

  z_aux(1) = zfield

  IF(l_fast_sum) THEN
    global_sum = simple_sum(z_aux, SIZE(z_aux), p_comm_glob, opt_iroot)
  ELSE
    ! Note: For (l_fast_sum == .FALSE.) there is no special
    !       implementation, which gathers only at rank
    !       "iroot". Instead, we always do an ALLREDUCE here.
    global_sum = order_insensit_ieee64_sum(z_aux, SIZE(z_aux), p_comm_glob)
  ENDIF

  IF(p_test_run .AND. do_sync_checks) THEN
    IF(l_fast_sum) THEN
      CALL check_result( (/ global_sum /), 'global_sum_array', sum_on_testpe)
      global_sum = sum_on_testpe(1)
    ELSE
      CALL check_result( (/ global_sum /), 'global_sum_array')
    ENDIF
  ENDIF
END FUNCTION global_sum_array_0d

!-------------------------------------------------------------------------
!! Calculates the global sum of zfield and checks for consistency
!! when doing a verification run.
!! This routine should be called outside an OMP parallel Region!
!!
FUNCTION global_sum_array_1d (zfield, opt_iroot) RESULT (global_sum)

  REAL(wp),          INTENT(in) :: zfield(:)
  INTEGER,  INTENT(IN),OPTIONAL :: opt_iroot
  REAL(wp)                      :: global_sum
  REAL(wp)                      :: sum_on_testpe(1)

  INTEGER :: p_comm_glob
!-----------------------------------------------------------------------

  IF(comm_lev==0) THEN
    p_comm_glob = p_comm_work
  ELSE
    p_comm_glob = glob_comm(comm_lev)
  ENDIF

  IF(l_fast_sum) THEN
    global_sum = simple_sum(zfield, SIZE(zfield), p_comm_glob, opt_iroot)
  ELSE
    global_sum = order_insensit_ieee64_sum(zfield, SIZE(zfield), p_comm_glob)
  ENDIF

  IF(p_test_run .AND. do_sync_checks) THEN
    IF(l_fast_sum) THEN
      CALL check_result( (/ global_sum /), 'global_sum_array', sum_on_testpe)
      global_sum = sum_on_testpe(1)
    ELSE
      CALL check_result( (/ global_sum /), 'global_sum_array')
    ENDIF
  ENDIF


END FUNCTION global_sum_array_1d


!-------------------------------------------------------------------------
!! Calculates the global sum of zfield and checks for consistency
!! when doing a verification run.
!! This routine should be called outside an OMP parallel Region!
!!
FUNCTION global_sum_array_2d (zfield, lacc) RESULT (global_sum)

  REAL(wp),          INTENT(in) :: zfield(:, :)
  LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc
  REAL(wp)                      :: global_sum
  REAL(wp)                      :: sum_on_testpe(1)

  INTEGER :: p_comm_glob
!-----------------------------------------------------------------------

  IF(comm_lev==0) THEN
    p_comm_glob = p_comm_work
  ELSE
    p_comm_glob = glob_comm(comm_lev)
  ENDIF

  IF(l_fast_sum) THEN
    CALL assert_acc_host_only("mo_sync:global_sum_array_2d-l_fast_sum", lacc)
    global_sum = simple_sum(zfield, SIZE(zfield), p_comm_glob)
  ELSE
    global_sum = order_insensit_ieee64_sum(zfield, SIZE(zfield), p_comm_glob, lacc=lacc)
  ENDIF

  IF(p_test_run .AND. do_sync_checks) THEN
    CALL assert_acc_host_only("mo_sync:global_sum_array_2d-do_sync_checks", lacc)
    IF(l_fast_sum) THEN
      CALL check_result( (/ global_sum /), 'global_sum_array', sum_on_testpe)
      global_sum = sum_on_testpe(1)
    ELSE
      CALL check_result( (/ global_sum /), 'global_sum_array')
    ENDIF
  ENDIF


END FUNCTION global_sum_array_2d

!-------------------------------------------------------------------------
!! Calculates the global sum of zfield and checks for consistency
!! when doing a verification run.
!! This routine should be called outside an OMP parallel Region!
!!
FUNCTION global_sum_array_3d (zfield) RESULT (global_sum)

  REAL(wp),          INTENT(in) :: zfield(:,:,:)
  REAL(wp)                      :: global_sum
  REAL(wp)                      :: sum_on_testpe(1)

  INTEGER :: p_comm_glob
!-----------------------------------------------------------------------

  IF(comm_lev==0) THEN
    p_comm_glob = p_comm_work
  ELSE
    p_comm_glob = glob_comm(comm_lev)
  ENDIF

  IF(l_fast_sum) THEN
    global_sum = simple_sum(zfield, SIZE(zfield), p_comm_glob)
  ELSE
    global_sum = order_insensit_ieee64_sum(zfield, SIZE(zfield), p_comm_glob)
  ENDIF

  IF(p_test_run .AND. do_sync_checks) THEN
    IF(l_fast_sum) THEN
      CALL check_result( (/ global_sum /), 'global_sum_array', sum_on_testpe)
      global_sum = sum_on_testpe(1)
    ELSE
      CALL check_result( (/ global_sum /), 'global_sum_array')
    ENDIF
  ENDIF

END FUNCTION global_sum_array_3d
!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
!! Calculates the global sum of zfield and checks for consistency
!! when doing a verification run.
!! This routine should be called from within an OMP parallel Region!
!!
FUNCTION omp_global_sum_array_1d (zfield) RESULT (global_sum)

!
  REAL(wp),          INTENT(in) :: zfield(:)
  REAL(wp)                      :: global_sum
  REAL(wp)                      :: sum_on_testpe(1)

  INTEGER :: p_comm_glob
!-----------------------------------------------------------------------

  IF(comm_lev==0) THEN
    p_comm_glob = p_comm_work
  ELSE
    p_comm_glob = glob_comm(comm_lev)
  ENDIF

  IF(l_fast_sum) THEN
    global_sum = omp_simple_sum(zfield, SIZE(zfield), p_comm_glob)
  ELSE
    global_sum = omp_order_insensit_ieee64_sum(zfield, SIZE(zfield), p_comm_glob)
  ENDIF

  IF(p_test_run .AND. do_sync_checks) THEN
!$OMP BARRIER
!$OMP MASTER
    IF(l_fast_sum) THEN
      CALL check_result( (/ global_sum /), 'global_sum_array', sum_on_testpe)
      global_sum = sum_on_testpe(1)
    ELSE
      CALL check_result( (/ global_sum /), 'global_sum_array')
    ENDIF
!$OMP END MASTER
!$OMP BARRIER
  ENDIF

END FUNCTION omp_global_sum_array_1d

!-------------------------------------------------------------------------
!! Calculates the global sum of zfield and checks for consistency
!! when doing a verification run.
!! This routine should be called from within an OMP parallel Region!
!!
FUNCTION omp_global_sum_array_2d (zfield) RESULT (global_sum)

!
  REAL(wp),          INTENT(in) :: zfield(:, :)
  REAL(wp)                      :: global_sum
  REAL(wp)                      :: sum_on_testpe(1)

  INTEGER :: p_comm_glob
!-----------------------------------------------------------------------

  IF(comm_lev==0) THEN
    p_comm_glob = p_comm_work
  ELSE
    p_comm_glob = glob_comm(comm_lev)
  ENDIF

  IF(l_fast_sum) THEN
    global_sum = omp_simple_sum(zfield, SIZE(zfield), p_comm_glob)
  ELSE
    global_sum = omp_order_insensit_ieee64_sum(zfield, SIZE(zfield), p_comm_glob)
  ENDIF

  IF(p_test_run .AND. do_sync_checks) THEN
!$OMP BARRIER
!$OMP MASTER
    IF(l_fast_sum) THEN
      CALL check_result( (/ global_sum /), 'global_sum_array', sum_on_testpe)
      global_sum = sum_on_testpe(1)
    ELSE
      CALL check_result( (/ global_sum /), 'global_sum_array')
    ENDIF
!$OMP END MASTER
!$OMP BARRIER
  ENDIF

END FUNCTION omp_global_sum_array_2d


!-------------------------------------------------------------------------
!
!

!! Calculates the global sum of zfield and checks for consistency.
!!
!! Calculates the global sum of zfield and checks for consistency
!! when doing a verification run.
!! This routine has to be called from outside an OMP parallel Region!
!!
FUNCTION global_sum_array2 (zfield) RESULT (global_sum)

!
   REAL(wp), INTENT(in) :: zfield(:,:)
   REAL(wp)             :: global_sum

   INTEGER(i8)       :: itmp(2), isum(2), ival1, ival2
   INTEGER(i8)       :: ilsum(SIZE(zfield,2),2)
   INTEGER           :: i, j, iexp, nblks, nlen
   REAL(dp)          :: abs_max
   REAL(dp)          :: fact, r_fact, rval
   REAL(dp)          :: max_aux(SIZE(zfield,2))

   REAL(dp), PARAMETER :: two_40 = 1099511627776._dp ! 2.**40
   REAL(dp), PARAMETER :: r_two_40 = 1._dp/two_40

   INTEGER(i8), PARAMETER :: mask40 = INT(z'000000ffffffffff',i8)
   INTEGER :: p_comm_glob

!-----------------------------------------------------------------------

   ilsum = 0_i8
   nblks = SIZE(zfield,2)
   nlen  = SIZE(zfield,1)

   IF(comm_lev==0) THEN
     p_comm_glob = p_comm_work
   ELSE
     p_comm_glob = glob_comm(comm_lev)
   ENDIF

   ! Get the maximum absolute value of all numbers.
!$OMP PARALLEL PRIVATE(rval)
!$OMP DO PRIVATE(j)
   DO j=1,nblks
      max_aux(j) = MAXVAL(ABS(zfield(:,j)))
   ENDDO
!$OMP END DO

!$OMP MASTER
   rval = MAXVAL(max_aux)
   abs_max = p_max(rval, comm=p_comm_glob)

   ! Get the exponent of abs_max for scaling
   iexp = EXPONENT(abs_max)

   ! Calculate a factor for scaling the input numbers
   ! so that the maximum absolute value of a scaled number
   ! is below 2**40

   fact = SCALE(1._dp,40-iexp) ! same as 2**(40-iexp)
   r_fact = SCALE(1._dp,iexp-40) ! 1./fact

!$OMP END MASTER
!$OMP BARRIER

   ! Sum up all numbers as scaled integers

!$OMP DO PRIVATE(i, j, rval, ival1, ival2)
   DO j=1,nblks
     DO i=1,nlen
       ! Scale number into range -2**40 < rval < 2**40
       ! and store integer part in ival1

       rval = zfield(i,j)*fact
       ival1 = INT(rval,i8)

       ! Scale fraction by 2**40 and store integer part in ival2

       ival2 = INT((rval - REAL(ival1,dp))*two_40,i8)

       ! Sum up ival1 and ival2; since we are using 8-byte integers
       ! for the sum there are no problems with overflow

       ilsum(j,1) = ilsum(j,1) + ival1
       ilsum(j,2) = ilsum(j,2) + ival2
     ENDDO
   ENDDO
!$OMP END DO
!$OMP END PARALLEL

   ! If the exponent is too small, return 0 in order to avoid
   ! problems with overflow below

   IF(iexp < -980) THEN
      global_sum = 0._dp
      RETURN
   ENDIF

   itmp(1) = SUM(ilsum(:,1))
   itmp(2) = SUM(ilsum(:,2))
   isum = p_sum(itmp, comm=p_comm_glob)

   ! Scale integer numbers back to real numbers and add them.
   ! For safety, we use only positive INTEGERS < 2**40 when converting to REAL

    IF (isum(1) >= 0_i8) THEN
       ival1 = ISHFT(isum(1),-40)
       ival2 = IAND (isum(1),mask40)
       global_sum = (REAL(ival1,dp)*r_fact)*two_40 + REAL(ival2,dp)*r_fact
    ELSE
       ival1 = ISHFT(ABS(isum(1)),-40)
       ival2 = IAND (ABS(isum(1)),mask40)
       global_sum = (REAL(ival1,dp)*r_fact)*two_40 + REAL(ival2,dp)*r_fact
       global_sum = -global_sum
    ENDIF

    IF (isum(2) >= 0_i8) THEN
       ival1 = ISHFT(isum(2),-40)
       ival2 = IAND (isum(2),mask40)
       global_sum = global_sum + (REAL(ival1,dp)*r_fact) + (REAL(ival2,dp)*r_fact)*r_two_40
    ELSE
       ival1 = ISHFT(ABS(isum(2)),-40)
       ival2 = IAND (ABS(isum(2)),mask40)
       global_sum = global_sum - (REAL(ival1,dp)*r_fact) - (REAL(ival2,dp)*r_fact)*r_two_40
    ENDIF

   IF(p_test_run .AND. do_sync_checks) CALL check_result( (/ global_sum /), 'global_sum_array2')

END FUNCTION global_sum_array2


!! Calculates the global sum of 3D and/or 4D input fields and checks for
!! consistency when doing a verification run.
!! This routine has to be called from outside an OMP parallel Region!
!! If ldiff=.true., a difference or ratio between f??in and f??d is computed
!! nfields specifies the number of 3D input fields, or pairs of input fields if
!! ldiff=.true. (when providing a 4D field as input, it counts as "n" 3D fields)
!! diffmask specifies if differences or ratios are to be computed.
!!
FUNCTION global_sum_array3 (nfields,ldiff,f3din,f3dd,f3din2,f3dd2,f4din,f4dd,diffmask) &
         RESULT (global_sum)

   REAL(wp), INTENT(in), TARGET           :: f3din(:,:,:)
   REAL(wp), INTENT(in), TARGET, OPTIONAL :: f3dd(:,:,:), f3din2(:,:,:), f3dd2(:,:,:), &
                                             f4din(:,:,:,:), f4dd(:,:,:,:)

   INTEGER,  INTENT(in)           :: nfields  ! Total number of 3D input fields
   LOGICAL,  INTENT(in)           :: ldiff    ! .TRUE.: compute differences or ratios (depending
                                              ! on diffmask) between f??in and f??d
   INTEGER,  INTENT(in), OPTIONAL :: diffmask(nfields) ! 1: compute differences, 2: compute ratios

   REAL(wp)                       :: global_sum(nfields*SIZE(f3din,2))

   INTEGER(i8)       :: ival1, ival2
   INTEGER           :: i, j, k, kk, n, nblks, nlen, nlev, nblks2, nfldtot, &
                        nlevtot, icnt, icnt2, nnblks(2*nfields)
   REAL(dp)          :: rlog2, rval1

   INTEGER(i8), ALLOCATABLE :: itmp(:,:),isum(:,:)
   INTEGER    , ALLOCATABLE :: iexp(:)
   REAL(dp),    ALLOCATABLE :: fact(:),r_fact(:),rval(:),abs_max(:)
   REAL(wp),    ALLOCATABLE :: aux_sum(:)

   TYPE t_fieldptr
     REAL(wp), POINTER :: fld(:,:,:)
   END TYPE t_fieldptr

   TYPE(t_fieldptr), ALLOCATABLE :: ff(:)

   REAL(dp), PARAMETER :: two_40 = 1099511627776._dp ! 2.**40
   REAL(dp), PARAMETER :: r_two_40 = 1._dp/two_40

   INTEGER(i8), PARAMETER :: mask40 = INT(z'000000ffffffffff',i8)
   INTEGER :: p_comm_glob

!-----------------------------------------------------------------------

   IF (ldiff .AND. .NOT.(PRESENT(diffmask))) THEN
     CALL finish('global_sum_array3','ldiff=.TRUE. requires the presence of diffmask')
   ENDIF

   IF(comm_lev==0) THEN
     p_comm_glob = p_comm_work
   ELSE
     p_comm_glob = glob_comm(comm_lev)
   ENDIF

   rlog2 = LOG(2._dp)
   nlen  = SIZE(f3din,1)
   nlev  = SIZE(f3din,2)
   nblks = SIZE(f3din,3)
   IF (ldiff) THEN
     nblks2 = SIZE(f3dd,3)
     nfldtot = 2*nfields
   ELSE
     nblks2 = 0
     nfldtot = nfields
   ENDIF
   nlevtot           = nlev*nfldtot
   nnblks(1:nfields) = nblks

   ALLOCATE (ff(nfldtot),itmp(nlevtot,2),isum(nlevtot,2),iexp(nlevtot),fact(nlevtot),&
             r_fact(nlevtot),rval(nlevtot),abs_max(nlevtot),aux_sum(nlevtot))

   itmp = 0_i8

   ! Set pointers to input fields
   icnt = 1
   icnt2 = 1
   ff(icnt)%fld => f3din
   IF (PRESENT(f3din2)) THEN
     icnt = icnt + 1
     icnt2 = icnt2 + 1
     ff(icnt)%fld => f3din2
   ENDIF
   IF (PRESENT(f4din)) THEN
     DO i = 1, nfields-icnt2
       icnt = icnt + 1
       ff(icnt)%fld => f4din(:,:,:,i)
     ENDDO
   ENDIF

   IF (ldiff) THEN
     icnt = icnt + 1
     icnt2 = 1
     ff(icnt)%fld => f3dd
     nnblks(icnt) = nblks2
     IF (PRESENT(f3dd2)) THEN
       icnt = icnt + 1
       icnt2 = icnt2 + 1
       ff(icnt)%fld => f3dd2(:,:,:)
       nnblks(icnt) = nblks2
     ENDIF
     IF (PRESENT(f4dd)) THEN
       DO i = 1, nfields-icnt2
         icnt = icnt + 1
         ff(icnt)%fld => f4dd(:,:,:,i)
         nnblks(icnt) = nblks2
       ENDDO
     ENDIF
   ENDIF

! OpenMP parallelization is done over the vertical levels in this routine
!$OMP PARALLEL
!$OMP DO PRIVATE(k,kk,n)
   ! Get the maximum absolute value of all field elements of a model level.
   DO k=1,nlevtot
     n = INT((k-1)/nlev)+1
     kk = k-(n-1)*nlev
     rval(k) = MAXVAL(ABS(ff(n)%fld(:,kk,:)))
   ENDDO
!$OMP END DO
!$OMP END PARALLEL

   abs_max = p_max(rval, comm=p_comm_glob)

   DO k=1,nlevtot
     abs_max(k) = MAX(abs_max(k),1.d-250)
     ! Get the exponent of abs_max for scaling
     iexp(k) = INT(LOG(abs_max(k))/rlog2)
     ! Calculate a factor for scaling the input numbers
     ! so that the maximum absolute value of a scaled number
     ! is below 2**40
     fact(k) = 2._wp**(40._wp-REAL(iexp(k),wp))
     r_fact(k) = 1._dp/fact(k)
   ENDDO

!$OMP PARALLEL
!$OMP DO PRIVATE(i, j, k, kk, n, ival1, rval1)
   DO k=1,nlevtot
     n = INT((k-1)/nlev)+1
     kk = k-(n-1)*nlev
     ! Sum up all numbers as scaled integers
!CDIR UNROLL=4
     DO j=1,nnblks(n)
       DO i=1,nlen
       ! Scale number into range -2**40 < rval < 2**40
       ! and store integer part in ival1

         rval1 = ff(n)%fld(i,kk,j)*fact(k)
         ival1 = INT(rval1,i8)

         ! Sum up ival1 and ival2; since we are using 8-byte integers
         ! for the sum there are no problems with overflow

         itmp(k,1) = itmp(k,1) + ival1
         itmp(k,2) = itmp(k,2) + INT((rval1 - REAL(ival1,dp))*two_40,i8)
       ENDDO
     ENDDO

   ENDDO
!$OMP END DO
!$OMP END PARALLEL

   isum(:,1) = p_sum(itmp(:,1), comm=p_comm_glob)
   isum(:,2) = p_sum(itmp(:,2), comm=p_comm_glob)

   DO k=1,nlevtot

     ! If the exponent is too small, return 0 in order to avoid
     ! problems with overflow below
     IF(iexp(k) < -980) THEN
       aux_sum(k) = 0._dp
       CYCLE
     ENDIF

     ! Scale integer numbers back to real numbers and add them.
     ! For safety, we use only positive INTEGERS < 2**40 when converting to REAL

     IF (isum(k,1) >= 0_i8) THEN
       ival1 = ISHFT(isum(k,1),-40)
       ival2 = IAND (isum(k,1),mask40)
       aux_sum(k) = (REAL(ival1,dp)*r_fact(k))*two_40 + REAL(ival2,dp)*r_fact(k)
     ELSE
       ival1 = ISHFT(ABS(isum(k,1)),-40)
       ival2 = IAND (ABS(isum(k,1)),mask40)
       aux_sum(k) = (REAL(ival1,dp)*r_fact(k))*two_40 + REAL(ival2,dp)*r_fact(k)
       aux_sum(k) = -aux_sum(k)
     ENDIF

     IF (isum(k,2) >= 0_i8) THEN
       ival1 = ISHFT(isum(k,2),-40)
       ival2 = IAND (isum(k,2),mask40)
       aux_sum(k) = aux_sum(k) + (REAL(ival1,dp)*r_fact(k)) + &
                    (REAL(ival2,dp)*r_fact(k))*r_two_40
     ELSE
       ival1 = ISHFT(ABS(isum(k,2)),-40)
       ival2 = IAND (ABS(isum(k,2)),mask40)
       aux_sum(k) = aux_sum(k) - (REAL(ival1,dp)*r_fact(k)) - &
                    (REAL(ival2,dp)*r_fact(k))*r_two_40
     ENDIF

   ENDDO

   IF (ldiff) THEN ! Compute difference or ratio of global sums
     n = nlevtot/2
     DO k=1,n
       i = INT((k-1)/nlev)+1
       IF (diffmask(i) == 1) THEN
         global_sum(k) = aux_sum(k) - aux_sum(k+n)
       ELSE
         global_sum(k) = MIN(1.001_wp, aux_sum(k) / MAX(1.e-50_wp,aux_sum(k+n)))
         global_sum(k) = MAX(0.999_wp, global_sum(k))
       ENDIF
     ENDDO
   ELSE
     n = nlevtot
     global_sum(1:n) = aux_sum(1:n)
   ENDIF

   IF(p_test_run .AND. do_sync_checks) CALL check_result(global_sum, 'global_sum_array3')

   DEALLOCATE (ff,itmp,isum,iexp,fact,r_fact,rval,abs_max,aux_sum)

END FUNCTION global_sum_array3
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> global_min/global_max family:
!! Calculate global min/max (using current communicator) and compare
!! result to result on test PE (in case of a test run)
FUNCTION global_min_0d(zfield) RESULT(global_min)

  REAL(wp), INTENT(IN) :: zfield
  REAL(wp) :: global_min

  IF(comm_lev==0) THEN
    global_min = p_min(zfield, comm=p_comm_work)
  ELSE
    global_min = p_min(zfield, comm=glob_comm(comm_lev))
  ENDIF

  IF(p_test_run .AND. do_sync_checks) CALL check_result( (/ global_min /), 'global_min' )

END FUNCTION global_min_0d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> global_min_0d for INTEGER
FUNCTION global_min_0di(zfield) RESULT(global_min)

  INTEGER, INTENT(IN) :: zfield
  INTEGER  :: global_min
  REAL(wp) :: global_min_check

  IF(comm_lev==0) THEN
    global_min = p_min(zfield, comm=p_comm_work)
  ELSE
    global_min = p_min(zfield, comm=glob_comm(comm_lev))
  ENDIF

  IF(p_test_run .AND. do_sync_checks) THEN
    global_min_check = REAL(global_min)
    CALL check_result( (/ global_min_check /), 'global_min' )
  ENDIF

END FUNCTION global_min_0di
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
FUNCTION global_min_1d(zfield) RESULT(global_min)

  REAL(wp), INTENT(IN) :: zfield(:)
  REAL(wp) :: global_min(SIZE(zfield))

  IF(comm_lev==0) THEN
    global_min = p_min(zfield, comm=p_comm_work)
  ELSE
    global_min = p_min(zfield, comm=glob_comm(comm_lev))
  ENDIF

  IF(p_test_run .AND. do_sync_checks) CALL check_result( global_min, 'global_min' )

END FUNCTION global_min_1d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!
! @param[out]   proc_id  (Optional:) PE number of maximum value
! @param[inout] keyval   (Optional:) additional meta information
! @param[in]    iroot    (Optional:) root PE, otherwise we perform an
!                                    ALL-TO-ALL operation
!
! The parameter @p keyval can be used to communicate
! additional data on the maximum value, e.g., the level
! index where the maximum occurred.
!
FUNCTION global_max_0d(zfield, proc_id, keyval, iroot, icomm) RESULT(global_max)

  REAL(wp), INTENT(IN) :: zfield
  INTEGER, OPTIONAL, INTENT(inout) :: proc_id
  INTEGER, OPTIONAL, INTENT(inout) :: keyval
  INTEGER, OPTIONAL, INTENT(in)    :: iroot                !< rank of collecting PE
  INTEGER, OPTIONAL, INTENT(in)    :: icomm                !< MPI communicator
  REAL(wp) :: global_max
  INTEGER  :: pcomm

  IF(comm_lev==0) THEN
    pcomm=p_comm_work
  ELSE
    pcomm=glob_comm(comm_lev)
  END IF
  IF (PRESENT(icomm)) pcomm = icomm

  IF (p_test_run) THEN ! all-to-all communication required
    global_max = p_max(zfield, proc_id=proc_id, keyval=keyval, comm=pcomm)
  ELSE
    global_max = p_max(zfield, proc_id=proc_id, keyval=keyval, comm=pcomm, root=iroot)
  ENDIF

  IF(p_test_run .AND. do_sync_checks) CALL check_result( (/ global_max /), 'global_max' )

END FUNCTION global_max_0d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> global_max_0d for INTEGER
FUNCTION global_max_0di(zfield) RESULT(global_max)

  INTEGER, INTENT(IN) :: zfield
  INTEGER  :: global_max
  REAL(wp) :: global_max_check
  INTEGER  :: pcomm

  IF(comm_lev==0) THEN
    pcomm=p_comm_work
  ELSE
    pcomm=glob_comm(comm_lev)
  END IF

  IF (p_test_run) THEN ! all-to-all communication required
    global_max = p_max(zfield, pcomm)
  ELSE
    global_max = p_max(zfield, pcomm)
  ENDIF

  IF(p_test_run .AND. do_sync_checks) THEN
    global_max_check = REAL(global_max)
    CALL check_result( (/ global_max_check /), 'global_max' )
  ENDIF

END FUNCTION global_max_0di
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!
! @param[out]   proc_id  (Optional:) PE number of maximum value
! @param[inout] keyval   (Optional:) additional meta information
! @param[in]    iroot    (Optional:) root PE, otherwise we perform an
!                                    ALL-TO-ALL operation
!
! The parameter @p keyval can be used to communicate
! additional data on the maximum value, e.g., the level
! index where the maximum occurred.
!
FUNCTION global_max_1d(zfield, proc_id, keyval, iroot, icomm) RESULT(global_max)

  REAL(wp), INTENT(IN) :: zfield(:)
  INTEGER, OPTIONAL, INTENT(inout) :: proc_id(SIZE(zfield))
  INTEGER, OPTIONAL, INTENT(inout) :: keyval(SIZE(zfield))
  INTEGER, OPTIONAL, INTENT(in)    :: iroot                !< rank of collecting PE
  INTEGER, OPTIONAL, INTENT(in)    :: icomm                !< MPI communicator
  REAL(wp) :: global_max(SIZE(zfield))
  INTEGER  :: pcomm

  IF(comm_lev==0) THEN
    pcomm=p_comm_work
  ELSE
    pcomm=glob_comm(comm_lev)
  END IF
  IF (PRESENT(icomm)) pcomm = icomm

  IF (p_test_run) THEN ! all-to-all communication required
    global_max = p_max(zfield, proc_id=proc_id, keyval=keyval, comm=pcomm)
  ELSE
    global_max = p_max(zfield, proc_id=proc_id, keyval=keyval, comm=pcomm, root=iroot)
  ENDIF

  IF(p_test_run .AND. do_sync_checks) CALL check_result( global_max, 'global_max' )

END FUNCTION global_max_1d

!-------------------------------------------------------------------------------

!> Checks if res is identical on Test PE and working PEs.
!! If res_on_testpe is present, it contains the result on the test PE on exit,
!! otherwise the routine finishes if the results are not identical.

SUBROUTINE check_result(res, routine, res_on_testpe)

  REAL(wp), INTENT(IN) :: res(:)
  CHARACTER(len=*), INTENT(IN) :: routine
  REAL(wp), INTENT(out), OPTIONAL :: res_on_testpe(:)

  REAL(wp) :: aux(SIZE(res))
  INTEGER :: k
  LOGICAL :: out_of_sync, is_mpi_test

  is_mpi_test = my_process_is_mpi_test()
  IF (is_mpi_test) THEN
    aux(:) = res(:)
  ELSE
    aux(:) = 0.0_wp ! Safety only
  END IF

  IF(comm_lev==0) THEN
    CALL p_bcast(aux, process_mpi_all_test_id, comm=p_comm_work_test)
  ELSE
    IF(get_my_mpi_all_id() == process_mpi_all_test_id) THEN
      CALL p_send(aux, comm_proc0(comm_lev)+p_work_pe0, 1)
    ELSE
      IF(p_pe_work==comm_proc0(comm_lev)) CALL p_recv(aux, process_mpi_all_test_id, 1)
      CALL p_bcast(aux, 0, comm=glob_comm(comm_lev))
    ENDIF
  ENDIF

  IF( .NOT. is_mpi_test .AND. l_log_checks .AND. log_unit>0) THEN
    DO k = 1, SIZE(res)
      WRITE(log_unit,'(a,2g25.18,a,g25.18)') routine,aux(k),res(k),&
           ' Error: ',ABS(aux(k)-res(k))
    END DO
  END IF
  IF(PRESENT(res_on_testpe)) THEN
    res_on_testpe = aux
  ELSE
    out_of_sync = .FALSE.
    DO k = 1, SIZE(res)
      ! Check if result is identical
      out_of_sync = out_of_sync .OR. (aux(k)/=res(k))
    END DO
    IF (out_of_sync) CALL finish(routine, 'Result out of sync')
  END IF


END SUBROUTINE check_result


!-------------------------------------------------------------------------------
!! This routine calculates the sum of an array of IEEE 64 bit
!! floating point values in an order insensitve way.
!! This is done by calculating the sum in INTEGER arithmetic,
!! so it is always exactly the same for any permutation of the numbers
!! (even among several processors).
!! Since 60 bits are used for the mantissa of each operand (and even more
!! for the accumulator), this sum should be almost always more precise
!! than naivly summing up the operands.
!! ATTENTION: When compiled with OpenMP in effect, this routine
!! should be called from a parallel region!!!!
!!
FUNCTION omp_order_insensit_ieee64_sum(vals, num_vals, mpi_comm) RESULT(global_sum)

!
   INTEGER  :: num_vals, mpi_comm
   REAL(dp) :: vals(num_vals)

   REAL(dp) :: global_sum

   INTEGER(i8)       :: itmp(2),  ival1, ival2
   INTEGER(i8), SAVE :: isum(2) ! This must be a shared variable
   INTEGER           :: i, iexp
   REAL(dp), SAVE    :: abs_max ! This must be a shared variable
   REAL(dp)          :: fact, r_fact, rval

   REAL(dp), PARAMETER :: two_30 = 1073741824._dp ! 2.**30
   REAL(dp), PARAMETER :: r_two_30 = 1._dp/two_30

   INTEGER(i8), PARAMETER :: mask30 = INT(z'000000003fffffff',i8)

!-----------------------------------------------------------------------
!    IF (activate_sync_timers) CALL timer_start(timer_omp_ordglb_sum)

   ! Set shared variables in a MASTER region

!$OMP MASTER
   abs_max = 0._dp
   isum(:) = 0_i8
!$OMP END MASTER
!$OMP BARRIER

   ! Get the maximum absolute value of all numbers.

!$OMP DO PRIVATE(i), REDUCTION(MAX:abs_max)
   DO i=1,num_vals
      abs_max = MAX(abs_max, ABS(vals(i)))
   ENDDO
!$OMP END DO

!$OMP MASTER
   rval = abs_max
   abs_max = p_max(rval, comm=mpi_comm)
!$OMP END MASTER
!$OMP BARRIER

   ! If abs_max is 0, all input values are 0
   ! and we are done

   IF(abs_max == 0.0_dp) THEN
      global_sum = 0._dp
      RETURN
   ENDIF

   ! Get the exponent of abs_max for scaling

   iexp = EXPONENT(abs_max)

   ! If the exponent is too small, return 0 in order to avoid
   ! problems with overflow below

   IF(iexp < -980) THEN
      global_sum = 0._dp
      RETURN
   ENDIF

   ! Calculate a factor for scaling the input numbers
   ! so that the maximum absolute value of a scaled number
   ! is below 2**30

   fact = SCALE(1._dp,30-iexp) ! same as 2**(30-iexp)
   r_fact = SCALE(1._dp,iexp-30) ! 1./fact

   ! Sum up all numbers as scaled integers

!$OMP DO PRIVATE(i, rval, ival1, ival2), REDUCTION(+:isum)
   DO i=1,num_vals

      ! Scale number into range -2**30 < rval < 2**30
      ! and store integer part in ival1

      rval = vals(i)*fact
      ival1 = INT(rval,i8)

      ! Scale fraction by 2**30 and store integer part in ival2

      ival2 = INT((rval - REAL(ival1,dp))*two_30,i8)

      ! Sum up ival1 and ival2; since we are using 8-byte integers
      ! for the sum there are no problems with overflow

      isum(1) = isum(1) + ival1
      isum(2) = isum(2) + ival2

   ENDDO
!$OMP END DO

!$OMP MASTER
   itmp = isum
   isum = p_sum(itmp, comm=mpi_comm)
!$OMP END MASTER
!$OMP BARRIER

   ! Scale integer numbers back to real numbers and add them.
   ! For safety, we use only positive INTEGERS < 2**30 when converting to REAL

    IF(isum(1) >= 0_i8)THEN
       ival1 = ISHFT(isum(1),-30)
       ival2 = IAND (isum(1),mask30)
       global_sum = (REAL(ival1,dp)*r_fact)*two_30 + REAL(ival2,dp)*r_fact
    ELSE
       ival1 = ISHFT(ABS(isum(1)),-30)
       ival2 = IAND (ABS(isum(1)),mask30)
       global_sum = (REAL(ival1,dp)*r_fact)*two_30 + REAL(ival2,dp)*r_fact
       global_sum = -global_sum
    ENDIF

    IF(isum(2) >= 0_i8)THEN
       ival1 = ISHFT(isum(2),-30)
       ival2 = IAND (isum(2),mask30)
       global_sum = global_sum + (REAL(ival1,dp)*r_fact) + (REAL(ival2,dp)*r_fact)*r_two_30
    ELSE
       ival1 = ISHFT(ABS(isum(2)),-30)
       ival2 = IAND (ABS(isum(2)),mask30)
       global_sum = global_sum - (REAL(ival1,dp)*r_fact) - (REAL(ival2,dp)*r_fact)*r_two_30
    ENDIF

!    IF (activate_sync_timers) CALL timer_stop(timer_omp_ordglb_sum)

END FUNCTION omp_order_insensit_ieee64_sum
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!! This routine calculates the sum of an array of IEEE 64 bit
!! floating point values in an order insensitve way.
!! This is done by calculating the sum in INTEGER arithmetic,
!! so it is always exactly the same for any permutation of the numbers
!! (even among several processors).
!! Since 60 bits are used for the mantissa of each operand (and even more
!! for the accumulator), this sum should be almost always more precise
!! than naivly summing up the operands.
!! ATTENTION: When compiled with OpenMP in effect, this routine
!! should be called outside an omp parallel region!!!!
!! When used with OpenACC, the result is provided on host only.
!!
FUNCTION order_insensit_ieee64_sum(vals, num_vals, mpi_comm, lacc) RESULT(global_sum)

!
   INTEGER, INTENT(IN)  :: num_vals, mpi_comm
   REAL(dp), INTENT(IN) :: vals(num_vals)
   LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

   REAL(dp) :: global_sum

   INTEGER(i8)       :: itmp(2),  ival1, ival2
   INTEGER(i8)       :: isum_1, isum_2, isum(2)
   INTEGER           :: i, iexp
   REAL(dp)          :: abs_max
   REAL(dp)          :: fact, r_fact, rval
   LOGICAL :: lzacc ! non-optional version of lacc

   REAL(dp), PARAMETER :: two_30 = 1073741824._dp ! 2.**30
   REAL(dp), PARAMETER :: r_two_30 = 1._dp/two_30

   INTEGER(i8), PARAMETER :: mask30 = INT(z'000000003fffffff',i8)

!-----------------------------------------------------------------------
   IF (activate_sync_timers) CALL timer_start(timer_ordglb_sum)

   CALL set_acc_host_or_device(lzacc, lacc)

   ! Set shared variables in a MASTER region
   abs_max = 0._dp
   isum_1 = 0_i8
   isum_2 = 0_i8
   ! Get the maximum absolute value of all numbers.
   !$ACC PARALLEL DEFAULT(PRESENT) REDUCTION(MAX: abs_max) ASYNC(1) IF(lzacc)
   !$ACC LOOP GANG VECTOR REDUCTION(MAX: abs_max)
   DO i=1,num_vals
      abs_max = MAX(abs_max, ABS(vals(i)))
   ENDDO
   !$ACC END PARALLEL
   !$ACC WAIT IF(lzacc)
   rval = abs_max
   abs_max = p_max(rval, comm=mpi_comm)
   ! If abs_max is 0, all input values are 0
   ! and we are done
   IF(abs_max == 0.0_dp) THEN
      global_sum = 0._dp
      IF (activate_sync_timers) CALL timer_stop(timer_ordglb_sum)
      RETURN
   ENDIF

   ! Get the exponent of abs_max for scaling

   iexp = EXPONENT(abs_max)

   ! If the exponent is too small, return 0 in order to avoid
   ! problems with overflow below

   IF(iexp < -980) THEN
      global_sum = 0._dp
      IF (activate_sync_timers) CALL timer_stop(timer_ordglb_sum)
      RETURN
   ENDIF

   ! Calculate a factor for scaling the input numbers
   ! so that the maximum absolute value of a scaled number
   ! is below 2**30

   fact = SCALE(1._dp,30-iexp) ! same as 2**(30-iexp)
   r_fact = SCALE(1._dp,iexp-30) ! 1./fact

   ! Sum up all numbers as scaled integers
   !$ACC PARALLEL DEFAULT(PRESENT) REDUCTION(+: isum_1, isum_2) ASYNC(1) IF(lzacc)
   !$ACC LOOP GANG VECTOR PRIVATE(rval, ival1, ival2) REDUCTION(+: isum_1, isum_2)
   DO i=1,num_vals

      ! Scale number into range -2**30 < rval < 2**30
      ! and store integer part in ival1

      rval = vals(i)*fact
      ival1 = INT(rval,i8)

      ! Scale fraction by 2**30 and store integer part in ival2

      ival2 = INT((rval - REAL(ival1,dp))*two_30,i8)

      ! Sum up ival1 and ival2; since we are using 8-byte integers
      ! for the sum there are no problems with overflow

      isum_1 = isum_1 + ival1
      isum_2 = isum_2 + ival2

   ENDDO
   !$ACC END PARALLEL
   !$ACC WAIT IF(lzacc)

   itmp = (/ isum_1, isum_2 /)
   isum = p_sum(itmp, comm=mpi_comm)

   ! Scale integer numbers back to real numbers and add them.
   ! For safety, we use only positive INTEGERS < 2**30 when converting to REAL

    IF(isum(1) >= 0_i8)THEN
       ival1 = ISHFT(isum(1),-30)
       ival2 = IAND (isum(1),mask30)
       global_sum = (REAL(ival1,dp)*r_fact)*two_30 + REAL(ival2,dp)*r_fact
    ELSE
       ival1 = ISHFT(ABS(isum(1)),-30)
       ival2 = IAND (ABS(isum(1)),mask30)
       global_sum = (REAL(ival1,dp)*r_fact)*two_30 + REAL(ival2,dp)*r_fact
       global_sum = -global_sum
    ENDIF

    IF(isum(2) >= 0_i8)THEN
       ival1 = ISHFT(isum(2),-30)
       ival2 = IAND (isum(2),mask30)
       global_sum = global_sum + (REAL(ival1,dp)*r_fact) + (REAL(ival2,dp)*r_fact)*r_two_30
    ELSE
       ival1 = ISHFT(ABS(isum(2)),-30)
       ival2 = IAND (ABS(isum(2)),mask30)
       global_sum = global_sum - (REAL(ival1,dp)*r_fact) - (REAL(ival2,dp)*r_fact)*r_two_30
    ENDIF

    IF (activate_sync_timers) CALL timer_stop(timer_ordglb_sum)

END FUNCTION order_insensit_ieee64_sum
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!! This routine calculates the sum of an array in the
!! straightforward way in parallel.
!! ATTENTION: When compiled with OpenMP in effect, this routine
!! should be called outside a parallel omp region!!!!
!!
!! @param[in] opt_iroot (Optional:) root PE, otherwise we perform an
!!            ALL-TO-ALL operation.
!!
FUNCTION simple_sum(vals, num_vals, mpi_comm, opt_iroot) RESULT(global_sum)

  INTEGER                        :: num_vals, mpi_comm
  REAL(dp)                       :: vals(num_vals)
  INTEGER,  INTENT(IN), OPTIONAL :: opt_iroot
  REAL(dp)                       :: global_sum
  ! local variables
  INTEGER        :: i
  REAL(dp), SAVE :: s, res

!-----------------------------------------------------------------------
   IF (activate_sync_timers) CALL timer_start(timer_global_sum)

   s   = 0._dp
   res = 0._dp

   ! Sum up all numbers
   DO i=1,num_vals
      s = s + vals(i)
   ENDDO
   res = p_sum(s, comm=mpi_comm, root=opt_iroot)

   global_sum = res

   IF (activate_sync_timers) CALL timer_stop(timer_global_sum)

END FUNCTION simple_sum
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!! This routine calculates the sum of an array in the
!! straightforward way in parallel.
!! ATTENTION: When compiled with OpenMP in effect, this routine
!! should be called from a parallel region!!!!
!!
FUNCTION omp_simple_sum(vals, num_vals, mpi_comm) RESULT(global_sum)

!
   INTEGER  :: num_vals, mpi_comm
   REAL(dp) :: vals(num_vals)

   REAL(dp) :: global_sum

   INTEGER :: i
   REAL(dp), SAVE :: s, res

!-----------------------------------------------------------------------
  IF (activate_sync_timers) CALL timer_start(timer_omp_global_sum)

   ! Set shared variables in a MASTER region

!$OMP MASTER
   s = 0._dp
!$OMP END MASTER
!$OMP BARRIER

   ! Sum up all numbers

!$OMP DO PRIVATE(i), REDUCTION(+:s)
   DO i=1,num_vals
      s = s + vals(i)
   ENDDO
!$OMP END DO

!$OMP MASTER
   res = p_sum(s, comm=mpi_comm)
!$OMP END MASTER
!$OMP BARRIER

   global_sum = res

   IF (activate_sync_timers) CALL timer_stop(timer_omp_global_sum)

END FUNCTION omp_simple_sum
!-------------------------------------------------------------------------


!! Computes and prints summary information on the domain decomposition
!!
SUBROUTINE decomposition_statistics(p_patch)

   TYPE(t_patch), INTENT(INOUT) :: p_patch

   REAL(wp) :: cellstat(0:7),edgestat(0:6),vertstat(0:5), csmax(0:6),csmin(0:6),csavg(0:7), &
               esmax(0:6),esmin(0:6),esavg(0:6),vsmax(0:5),vsmin(0:5),vsavg(0:5),avglat,avglon
   INTEGER  :: i_nchdom, i, i_pe, max_nprecv, i1, i2, i1m, i2m
   INTEGER,  ALLOCATABLE :: nprecv_buf(:),displs(:),recvlist_buf(:)
   REAL(wp), ALLOCATABLE :: avglat_buf(:),avglon_buf(:)
   INTEGER :: np_recv
   INTEGER, ALLOCATABLE :: pelist_recv(:)

!-----------------------------------------------------------------------

   i_nchdom = MAX(1,p_patch%n_childdom)

   cellstat(0) = REAL(nproma*(p_patch%cells%end_blk(grf_bdywidth_c,1)-1) + &
                      p_patch%cells%end_idx(grf_bdywidth_c,1),wp)
   cellstat(1) = REAL(nproma*(p_patch%cells%end_blk(min_rlcell_int,i_nchdom)-1) + &
                      p_patch%cells%end_idx(min_rlcell_int,i_nchdom),wp) - cellstat(0)
   cellstat(2) = REAL(nproma*(p_patch%cells%end_blk(min_rlcell_int-1,i_nchdom)-1) + &
                      p_patch%cells%end_idx(min_rlcell_int-1,i_nchdom),wp) - cellstat(0)
   cellstat(3) = REAL(nproma*(p_patch%cells%end_blk(min_rlcell_int-2,i_nchdom)-1) + &
                      p_patch%cells%end_idx(min_rlcell_int-2,i_nchdom),wp) - cellstat(0)
   cellstat(4) = REAL(nproma*(p_patch%cells%end_blk(min_rlcell_int-2,i_nchdom)-1) + &
                      p_patch%cells%end_idx(min_rlcell_int-2,i_nchdom),wp)
   cellstat(5) = REAL(get_np_send(p_patch%comm_pat_c),wp)
   cellstat(6) = REAL(get_np_recv(p_patch%comm_pat_c),wp)

   ! The purpose of this is to compute average quantities only over those PEs
   ! that actually contain grid points of a given model domain (relevant in the case
   ! of processor splitting)
   IF (cellstat(1) > 0._wp) THEN
     cellstat(7) = 1._wp
   ELSE
     cellstat(7) = 0._wp
   ENDIF

   edgestat(0) = REAL(nproma*(p_patch%edges%end_blk(grf_bdywidth_e,1)-1) + &
                      p_patch%edges%end_idx(grf_bdywidth_e,1),wp)
   edgestat(1) = REAL(nproma*(p_patch%edges%end_blk(min_rledge_int,i_nchdom)-1) + &
                      p_patch%edges%end_idx(min_rledge_int,i_nchdom),wp) - edgestat(0)
   edgestat(2) = REAL(nproma*(p_patch%edges%end_blk(min_rledge_int-1,i_nchdom)-1) + &
                      p_patch%edges%end_idx(min_rledge_int-1,i_nchdom),wp) - edgestat(0)
   edgestat(3) = REAL(nproma*(p_patch%edges%end_blk(min_rledge_int-2,i_nchdom)-1) + &
                      p_patch%edges%end_idx(min_rledge_int-2,i_nchdom),wp) - edgestat(0)
   edgestat(4) = REAL(nproma*(p_patch%edges%end_blk(min_rledge_int-3,i_nchdom)-1) + &
                      p_patch%edges%end_idx(min_rledge_int-3,i_nchdom),wp) - edgestat(0)
   edgestat(5) = REAL(get_np_send(p_patch%comm_pat_e),wp)
   edgestat(6) = REAL(get_np_recv(p_patch%comm_pat_e),wp)

   vertstat(0) = REAL(nproma*(p_patch%verts%end_blk(grf_bdywidth_c,1)-1) + &
                      p_patch%verts%end_idx(grf_bdywidth_c,1),wp)
   vertstat(1) = REAL(nproma*(p_patch%verts%end_blk(min_rlvert_int,i_nchdom)-1) + &
                      p_patch%verts%end_idx(min_rlvert_int,i_nchdom),wp) - vertstat(0)
   vertstat(2) = REAL(nproma*(p_patch%verts%end_blk(min_rlvert_int-1,i_nchdom)-1) + &
                      p_patch%verts%end_idx(min_rlvert_int-1,i_nchdom),wp) - vertstat(0)
   vertstat(3) = REAL(nproma*(p_patch%verts%end_blk(min_rlvert_int-2,i_nchdom)-1) + &
                      p_patch%verts%end_idx(min_rlvert_int-2,i_nchdom),wp) - vertstat(0)
   vertstat(4) = REAL(get_np_send(p_patch%comm_pat_v),wp)
   vertstat(5) = REAL(get_np_recv(p_patch%comm_pat_v),wp)

   ! Question: how can I exclude PEs containing zero grid points of a model domain
   ! from global minimum computation?
   csmax = global_max(cellstat(0:6))
   csmin = global_min(cellstat(0:6))
   DO i = 0, 7
     csavg(i) = global_sum_array(cellstat(i))
   ENDDO
   csavg(0:6) = csavg(0:6)/MAX(1._wp,csavg(7))
   esmax = global_max(edgestat)
   esmin = global_min(edgestat)
   DO i = 0, 6
     esavg(i) = global_sum_array(edgestat(i))/MAX(1._wp,csavg(7))
   ENDDO
   vsmax = global_max(vertstat)
   vsmin = global_min(vertstat)
   DO i = 0, 5
     vsavg(i) = global_sum_array(vertstat(i))/MAX(1._wp,csavg(7))
   ENDDO

   WRITE(message_text,'(a,i4)') 'grid ',p_patch%id
     CALL message('Information on domain decomposition',TRIM(message_text))
   WRITE(message_text,'(i6)') NINT(csavg(7))
     CALL message('Number of compute PEs used for this grid',TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(csmax(0)),NINT(csmin(0)),csavg(0)
     CALL message('#   lateral boundary cells', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(csmax(1)),NINT(csmin(1)),csavg(1)
     CALL message('#         prognostic cells', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(csmax(2)),NINT(csmin(2)),csavg(2)
     CALL message('# cells up to halo level 1', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(csmax(3)),NINT(csmin(3)),csavg(3)
     CALL message('# cells up to halo level 2', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(csmax(4)),NINT(csmin(4)),csavg(4)
     CALL message('#              cells total', TRIM(message_text))

   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(esmax(0)),NINT(esmin(0)),esavg(0)
     CALL message('#   lateral boundary edges', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(esmax(1)),NINT(esmin(1)),esavg(1)
     CALL message('#         prognostic edges', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(esmax(2)),NINT(esmin(2)),esavg(2)
     CALL message('# edges up to halo level 1', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(esmax(3)),NINT(esmin(3)),esavg(3)
     CALL message('# edges up to halo level 2', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(esmax(4)),NINT(esmin(4)),esavg(4)
     CALL message('# edges up to halo level 3', TRIM(message_text))

   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(vsmax(0)),NINT(vsmin(0)),vsavg(0)
     CALL message('#   lateral boundary verts', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(vsmax(1)),NINT(vsmin(1)),vsavg(1)
     CALL message('#         prognostic verts', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(vsmax(2)),NINT(vsmin(2)),vsavg(2)
     CALL message('# verts up to halo level 1', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(vsmax(3)),NINT(vsmin(3)),vsavg(3)
     CALL message('# verts up to halo level 2', TRIM(message_text))
     CALL message('','')
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(csmax(5)),NINT(csmin(5)),csavg(5)
     CALL message('# send PEs (cells)', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(csmax(6)),NINT(csmin(6)),csavg(6)
     CALL message('# recv PEs (cells)', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(esmax(5)),NINT(esmin(5)),esavg(5)
     CALL message('# send PEs (edges)', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(esmax(6)),NINT(esmin(6)),esavg(6)
     CALL message('# recv PEs (edges)', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(vsmax(4)),NINT(vsmin(4)),vsavg(4)
     CALL message('# send PEs (verts)', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(vsmax(5)),NINT(vsmin(5)),vsavg(5)
     CALL message('# recv PEs (verts)', TRIM(message_text))
     CALL message('','')

   ! Stop if the current model domain lives on only one processor. MPI communication
   ! will not work in this case
   IF (.NOT. p_test_run .AND. NINT(csavg(7)) <= 1) &
     CALL finish('Bad use of processor splitting','This grid is processed by only one PE')


   IF (msg_level >= 14) THEN ! add list of cell recv PEs for each PE

     WRITE(message_text,'(a,i4)') 'grid ',p_patch%id
     CALL message('List of receive PEs (cells)', TRIM(message_text))
     max_nprecv = NINT(csmax(6))

     ! Compute average latitude and longitude over prognostic grid points (including nest boundary)
     avglat = 0._wp
     avglon = 0._wp
     i2m = p_patch%cells%end_blk(min_rlcell_int,i_nchdom)
     DO i2 = 1, i2m
       i1m = nproma
       if (i2 == i2m) i1m = p_patch%cells%end_idx(min_rlcell_int,i_nchdom)
       DO i1 = 1, i1m
         avglat = avglat + p_patch%cells%center(i1,i2)%lat
         avglon = avglon + p_patch%cells%center(i1,i2)%lon
       ENDDO
     ENDDO

     IF (cellstat(0)+cellstat(1) > 0._wp) THEN
       avglat = avglat/(cellstat(0)+cellstat(1))*180._wp/pi
       avglon = avglon/(cellstat(0)+cellstat(1))*180._wp/pi
     ELSE
       avglat = 0._wp
       avglon = 0._wp
     ENDIF

     ALLOCATE(nprecv_buf(p_n_work),displs(p_n_work),recvlist_buf(max_nprecv*p_n_work), &
              avglat_buf(p_n_work),avglon_buf(p_n_work), &
              pelist_recv(get_np_recv(p_patch%comm_pat_c)))

     np_recv = get_np_recv(p_patch%comm_pat_c)

     CALL p_gather(np_recv, nprecv_buf, 0, p_comm_work)

     displs(:) = (/ ( (i_pe-1)*max_nprecv, i_pe=1, p_n_work) /)

     CALL get_pelist_recv(p_patch%comm_pat_c, pelist_recv)

     CALL p_gatherv(pelist_recv, get_np_recv(p_patch%comm_pat_c), &
                    recvlist_buf, nprecv_buf, displs, 0, p_comm_work)

     CALL p_gather(avglat, avglat_buf, 0, p_comm_work)
     CALL p_gather(avglon, avglon_buf, 0, p_comm_work)

     IF (p_pe == 0) THEN
       WRITE(0,'(a)') 'PE, avg lat/lon, list of receive PEs'
       DO i_pe = 1, p_n_work
         i1 = (i_pe-1)*max_nprecv+1
         i2 = i1-1+MIN(18,nprecv_buf(i_pe))
         IF (nprecv_buf(i_pe) > 0) WRITE(0,'(i6,a,2f8.2,18i6)') i_pe-1,':', &
           avglat_buf(i_pe),avglon_buf(i_pe),recvlist_buf(i1:i2)
       ENDDO
     ENDIF

     DEALLOCATE(nprecv_buf,displs,recvlist_buf,avglat_buf,avglon_buf, pelist_recv)

   ENDIF

END SUBROUTINE decomposition_statistics

!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
!> Does boundary exchange for 3D cell-based fields, collecting as many
!  fields as possible before actually performing the sync.
!
SUBROUTINE cumulative_sync_patch_array(typ, p_patch, f3d)
  INTEGER,       INTENT(IN)            :: typ
  TYPE(t_patch), INTENT(IN),    TARGET :: p_patch
  REAL(wp),      INTENT(INOUT), TARGET :: f3d(:,:,:)
  ! local variables
  CHARACTER(*), PARAMETER :: routine = modname//"::cumulative_sync_patch_array"
  INTEGER :: idx

  ! add pointer to list of cumulative sync fields:
  idx = ncumul_sync(typ,p_patch%id) + 1
  IF (idx > MAX_CUMULATIVE_SYNC) CALL finish(routine, "Internal error!")
  ncumul_sync(typ, p_patch%id)            =  idx
  cumul_sync(typ, p_patch%id,idx)%f3d     => f3d
  cumul_sync(typ, p_patch%id,idx)%p_patch => p_patch
  IF (idx == MAX_CUMULATIVE_SYNC) THEN
!CDIR NOIEXPAND
    CALL complete_cumulative_sync(typ, p_patch%id)
  END IF
END SUBROUTINE cumulative_sync_patch_array


!-------------------------------------------------------------------------
!> If there are any pending "cumulative sync" operations: Complete them!
!
RECURSIVE SUBROUTINE complete_cumulative_sync(opt_typ, opt_patch_id)
  INTEGER, OPTIONAL, INTENT(IN) :: opt_typ
  INTEGER, OPTIONAL, INTENT(IN) :: opt_patch_id
  ! local variables
  CHARACTER(*), PARAMETER :: routine = modname//"::complete_cumulative_sync"
  TYPE(t_patch), POINTER :: p_patch
  INTEGER :: i

  ! first, handle the case that this routine has been called for all
  ! sync types or/and all patches:
  IF (.NOT. PRESENT(opt_typ)) THEN
    DO i=1,4
!CDIR NOIEXPAND
      CALL complete_cumulative_sync(i, opt_patch_id)
    END DO
    RETURN
  END IF
  IF (.NOT. PRESENT(opt_patch_id)) THEN
    DO i=1,max_dom
!CDIR NOIEXPAND
      CALL complete_cumulative_sync(opt_typ,i)
    END DO
    RETURN
  END IF
  ! now, call "sync_patch_array_mult"
  IF (PRESENT(opt_typ) .AND. (PRESENT(opt_patch_id))) THEN
    p_patch => cumul_sync(opt_typ,opt_patch_id,1)%p_patch
    SELECT CASE(ncumul_sync(opt_typ,opt_patch_id))
    CASE (0)
      RETURN
    CASE (1)
      CALL sync_patch_array_mult(opt_typ, p_patch, 1, cumul_sync(opt_typ,opt_patch_id,1)%f3d)
    CASE (2)
      CALL sync_patch_array_mult(opt_typ, p_patch, 2, cumul_sync(opt_typ,opt_patch_id, 1)%f3d, &
        &                        cumul_sync(opt_typ,opt_patch_id,2)%f3d)
    CASE (3)
      CALL sync_patch_array_mult(opt_typ, p_patch, 3, cumul_sync(opt_typ,opt_patch_id,1)%f3d,  &
        &                        cumul_sync(opt_typ,opt_patch_id,2)%f3d, cumul_sync(opt_typ,opt_patch_id,3)%f3d)
    CASE (4)
      CALL sync_patch_array_mult(opt_typ, p_patch, 4, cumul_sync(opt_typ,opt_patch_id,1)%f3d,  &
        &                        cumul_sync(opt_typ,opt_patch_id,2)%f3d, cumul_sync(opt_typ,opt_patch_id,3)%f3d, &
        &                        cumul_sync(opt_typ,opt_patch_id,4)%f3d)
    CASE (5)
      CALL sync_patch_array_mult(opt_typ, p_patch, 5, cumul_sync(opt_typ,opt_patch_id,1)%f3d,  &
        &                        cumul_sync(opt_typ,opt_patch_id,2)%f3d, cumul_sync(opt_typ,opt_patch_id,3)%f3d, &
        &                        cumul_sync(opt_typ,opt_patch_id,4)%f3d, cumul_sync(opt_typ,opt_patch_id,5)%f3d)
    CASE DEFAULT
      CALL finish(routine, "Internal error!")
    END SELECT
    ncumul_sync(opt_typ,opt_patch_id) = 0 ! reset sync counter
  END IF
END SUBROUTINE complete_cumulative_sync


END MODULE mo_sync

