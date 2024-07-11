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
MODULE mo_reader_sst_sic

  USE mo_kind,                    ONLY: dp, wp, i8
  USE mo_parallel_config,         ONLY: get_nproma
  USE mo_exception,               ONLY: finish
  USE mo_reader_abstract,         ONLY: t_abstract_reader
  USE mo_io_units,                ONLY: FILENAME_MAX
  USE mo_model_domain,            ONLY: t_patch
  USE mo_netcdf_errhandler,       ONLY: nf
  USE mo_netcdf
  USE mtime,                      ONLY: julianday, juliandelta, getjuliandayfromdatetime, &
       &                                datetime, newdatetime, deallocatedatetime,        &
       &                                OPERATOR(+), ASSIGNMENT(=),                       &
       &                                no_of_ms_in_a_day, no_of_ms_in_a_hour,            &
       &                                no_of_ms_in_a_minute, no_of_ms_in_a_second 
  USE mo_mpi,                     ONLY: my_process_is_stdio, my_process_is_mpi_workroot, &
       &                                process_mpi_root_id, p_comm_work, p_bcast
  USE mo_read_netcdf_distributed, ONLY: distrib_nf_open, distrib_read, distrib_nf_close, &
       &                                idx_lvl_blk
  USE mo_fortran_tools,           ONLY: t_ptr_3d
#ifdef _OPENACC
  USE mo_mpi,                     ONLY: i_am_accel_node
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_sst_sic_reader

  TYPE, EXTENDS(t_abstract_reader) :: t_sst_sic_reader

    TYPE(t_patch), POINTER      :: p_patch => NULL()
    CHARACTER(len=FILENAME_MAX) :: filename
    INTEGER                     :: fileid, dist_fileid
    LOGICAL                     :: lopened = .FALSE.

  CONTAINS

    PROCEDURE :: init            => sst_sic_init_reader
    PROCEDURE :: get_one_timelev => sst_sic_get_one_timelevel
    PROCEDURE :: get_times       => sst_sic_get_times
    PROCEDURE :: deinit          => sst_sic_deinit_reader

    PROCEDURE :: get_nblks       => sst_sic_get_nblks
    PROCEDURE :: get_npromz      => sst_sic_get_npromz

  END TYPE t_sst_sic_reader

  CHARACTER(len=*), PARAMETER :: modname = 'mo_reader_sst_sic'

CONTAINS

  SUBROUTINE sst_sic_init_reader(this, p_patch, filename)
    CLASS(t_sst_sic_reader),    INTENT(inout) :: this
    TYPE(t_patch),      TARGET, INTENT(in   ) :: p_patch
    CHARACTER(len=*),           INTENT(in   ) :: filename

    CHARACTER(len=*), PARAMETER :: routine = 'sst_sic_init_reader'

    this%filename = TRIM(filename)

    this%p_patch => p_patch

    IF (.NOT. this%lopened) THEN
      IF (my_process_is_mpi_workroot()) THEN
        CALL nf(nf90_open(this%filename, nf90_nowrite, this%fileid), routine)
      ENDIF
      this%dist_fileid = distrib_nf_open(TRIM(this%filename))
      this%lopened = .TRUE.
    ENDIF

  END SUBROUTINE sst_sic_init_reader

  SUBROUTINE sst_sic_get_times (this, times)

    CLASS(t_sst_sic_reader),      INTENT(inout) :: this
    TYPE(julianday), ALLOCATABLE, INTENT(  out) :: times(:)

    INTEGER                       :: tvid, tdid
    CHARACTER(len=NF90_MAX_NAME)    :: cf_timeaxis_string
    CHARACTER(len=:), ALLOCATABLE :: epoch
    CHARACTER(len=:), ALLOCATABLE :: base_timeaxis_unit
    TYPE(datetime), POINTER       :: epoch_datetime
    TYPE(julianday)               :: epoch_jd
    TYPE(juliandelta)             :: offset
    INTEGER(i8)                   :: time_multiplicator
    REAL(wp), ALLOCATABLE         :: times_read(:)
    INTEGER                       :: ntimes
    INTEGER                       :: i

    CHARACTER(len=*), PARAMETER :: routine = 'sst_sic_get_times'

    IF (my_process_is_mpi_workroot()) THEN

      CALL nf(nf90_inq_varid(this%fileid, "time", tvid), routine)
      CALL nf(nf90_inq_dimid(this%fileid, "time", tdid), routine)
      CALL nf(nf90_inquire_dimension(this%fileid, tdid, len = ntimes), routine)

      ALLOCATE(times_read(ntimes))

      CALL nf(nf90_get_var(this%fileid, tvid, times_read), routine)
      CALL nf(nf90_get_att(this%fileid,   tvid, "units",cf_timeaxis_string), routine)

    ENDIF

    CALL p_bcast(ntimes, process_mpi_root_id, p_comm_work)
    IF (.NOT. ALLOCATED(times_read)) THEN
      ALLOCATE(times_read(ntimes))
    ENDIF
    CALL p_bcast(times_read, process_mpi_root_id, p_comm_work)
    CALL p_bcast(cf_timeaxis_string, process_mpi_root_id, p_comm_work)

    CALL get_cf_timeaxis_desc(TRIM(cf_timeaxis_string), epoch, base_timeaxis_unit)

    epoch_datetime => newdatetime(epoch)
    CALL getJulianDayFromDatetime(epoch_datetime, epoch_jd)
    CALL deallocateDatetime(epoch_datetime)
    
    SELECT CASE (base_timeaxis_unit)
    CASE('days')
      time_multiplicator = no_of_ms_in_a_day
    CASE('hours')
      time_multiplicator = no_of_ms_in_a_hour
    CASE('minutes')
      time_multiplicator = no_of_ms_in_a_minute
    CASE('seconds')
      time_multiplicator = no_of_ms_in_a_second
    END SELECT
    
    ALLOCATE(times(ntimes))

    DO i = 1, ntimes
      offset%sign = '+'
      offset%day  = INT((time_multiplicator * times_read(i))/86400000.0_wp,i8)
      offset%ms   = NINT(MOD(time_multiplicator * times_read(i), 86400000.0_wp),i8)
      times(i) = epoch_jd + offset
    ENDDO

  END SUBROUTINE sst_sic_get_times

  SUBROUTINE sst_sic_get_one_timelevel(this, timelevel, varname, dat)
    CLASS(t_sst_sic_reader),   INTENT(inout) :: this
    INTEGER,               INTENT(in   ) :: timelevel
    CHARACTER(len=*),      INTENT(in   ) :: varname
    REAL(dp), ALLOCATABLE, INTENT(inout) :: dat(:,:,:,:)
    REAL(dp), ALLOCATABLE, TARGET        :: temp(:,:,:,:)
    TYPE(t_ptr_3d)                       :: tmp(1)
    ! We need to turn off OpenACC to do host-based MPI
    LOGICAL                              :: init_i_am_accel_node

    CHARACTER(len=*), PARAMETER :: routine = 'sst_sic_get_one_timelevel'
    

    ALLOCATE(temp(get_nproma(), 1, this%p_patch%nblks_c, 1))
    temp(:,:,:,:) = -1.0_dp
    IF (.NOT. this%lopened) THEN
      CALL finish(modname, '6 hourly SST/Seaice file not open!')
    END IF
    tmp(1)%p => temp(:,:,:,1)
#ifdef _OPENACC
    init_i_am_accel_node = i_am_accel_node
    i_am_accel_node = .FALSE.
#endif
    CALL distrib_read(this%dist_fileid, varname, tmp, &
         & (/this%p_patch%cells%dist_io_data/), edim=(/1/), dimo=idx_lvl_blk, &
         & start_ext_dim=(/timelevel/), end_ext_dim=(/timelevel/))
#ifdef _OPENACC
    i_am_accel_node = init_i_am_accel_node
#endif
    IF (ALLOCATED(dat)) THEN
      IF (.NOT. ALL( SHAPE(dat) .EQ. (/get_nproma(), 1, this%p_patch%nblks_c, 1/) )) THEN
        DEALLOCATE(dat)
      END IF
    END IF
    IF (.NOT. ALLOCATED(dat)) THEN
      CALL MOVE_ALLOC(temp, dat)
    ELSE
      dat(:,:,:,:) = temp(:,:,:,:)
    END IF

    CALL sst_sic_replace_missval(this, dat, -1.0_wp)

  END SUBROUTINE sst_sic_get_one_timelevel

  SUBROUTINE sst_sic_replace_missval (this, dat, new_missval)
    CLASS(t_sst_sic_reader), INTENT(inout) :: this
    REAL(wp),                INTENT(inout) :: dat(:,:,:,:)
    REAL(wp),                INTENT(in   ) :: new_missval

    WHERE (dat < -1e10_wp)
      dat = new_missval
    END WHERE

  END SUBROUTINE sst_sic_replace_missval

  FUNCTION sst_sic_get_nblks (this) RESULT(nblks)
    CLASS(t_sst_sic_reader), INTENT(in   ) :: this
    INTEGER                                :: nblks
    nblks = this%p_patch%nblks_c
  END FUNCTION sst_sic_get_nblks

  FUNCTION sst_sic_get_npromz (this) RESULT(npromz)
    CLASS(t_sst_sic_reader), INTENT(in   ) :: this
    INTEGER                                :: npromz
    npromz = this%p_patch%npromz_c
  END FUNCTION sst_sic_get_npromz

  SUBROUTINE sst_sic_deinit_reader(this)
    CLASS(t_sst_sic_reader), INTENT(inout) :: this
    CHARACTER(len=*), PARAMETER :: routine = 'sst_sic_deinit_reader'
    IF (ASSOCIATED(this%p_patch)) NULLIFY(this%p_patch)
    IF (this%lopened) THEN
      IF (my_process_is_mpi_workroot()) THEN
        CALL nf(nf90_close(this%fileid), routine)
      END IF
      CALL distrib_nf_close(this%dist_fileid)
    END IF
  END SUBROUTINE sst_sic_deinit_reader

  SUBROUTINE get_cf_timeaxis_desc(cf_timeaxis_string, epoch, base_timeaxis_unit)
    CHARACTER(len=*), INTENT(in) :: cf_timeaxis_string
    CHARACTER(len=:), ALLOCATABLE, INTENT(out) :: epoch, base_timeaxis_unit
    
    ! The CF convention allows for a timezone to be included. We will
    ! ignore that one for all , but gets stored to word(5), if
    ! provided, to keep the algorithm simple.

    CHARACTER(len=16) :: word(5)
    INTEGER :: pos1, pos2, n

    pos1 = 1; pos2 = 0; n = 0;
    word(:) = ""

    DO
      pos2 = INDEX(cf_timeaxis_string(pos1:), " ")
      IF (pos2 == 0) THEN
        n = n + 1
        word(n) = cf_timeaxis_string(pos1:)
        EXIT
      ENDIF
      n = n + 1
      word(n) = cf_timeaxis_string(pos1:pos1+pos2-2)
      pos1 = pos2+pos1
    ENDDO

 ! correct the date part
    normalize_date: BLOCK
      INTEGER :: idx1, idx2
      INTEGER :: year, month, day
      idx1 = INDEX(word(3), '-')    
      idx2 = INDEX(word(3)(idx1+1:), '-')+idx1
      READ(word(3)(      :idx1-1),*) year
      READ(word(3)(idx1+1:idx2-1),*) month
      READ(word(3)(idx2+1:      ),*) day
      WRITE(word(3),'(i0,a,i2.2,a,i2.2)') year, '-', month, '-', day
    END BLOCK normalize_date
    
    IF (word(4) /= "") THEN
      epoch = TRIM(word(3))//'T'//TRIM(word(4))
    ELSE
      epoch = TRIM(word(3))
    ENDIF

    base_timeaxis_unit = TRIM(word(1))

  END SUBROUTINE get_cf_timeaxis_desc
  
END MODULE mo_reader_sst_sic
