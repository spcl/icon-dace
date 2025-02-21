! This is  a clone of the respective ECHAM routine
!
! TODO: ctfreez in echam = 271.38, this is 271.45 K
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

MODULE mo_bc_sst_sic
  
  USE mo_kind,               ONLY: dp, i8
  USE mo_exception,          ONLY: finish, message, message_text
  USE mo_mpi,                ONLY: my_process_is_mpi_workroot, p_bcast, &
    &                              process_mpi_root_id, p_comm_work
  USE mo_model_domain,       ONLY: t_patch
  USE mo_grid_config,        ONLY: n_dom
  USE mo_parallel_config,    ONLY: nproma
  USE mo_physical_constants, ONLY: tf_salt !, tmelt
  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH, max_dom
  USE mo_cdi,                ONLY: streamOpenRead, streamInqVlist, streamClose, &
    & vlistInqTaxis, streamInqTimestep, taxisInqVdate, streamReadVarSlice
  USE mo_util_cdi,           ONLY: cdiGetStringError
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights, &
       &                               calculate_time_interpolation_weights
  USE mo_time_config,        ONLY: time_config
 
  IMPLICIT NONE

  PRIVATE

  TYPE t_ext_sea
    REAL(dp), CONTIGUOUS, POINTER :: sst(:,:,:) => NULL()
    REAL(dp), CONTIGUOUS, POINTER :: sic(:,:,:) => NULL()
  END TYPE t_ext_sea

  TYPE(t_ext_sea), TARGET :: ext_sea(max_dom)

  PUBLIC :: read_bc_sst_sic
  PUBLIC :: bc_sst_sic_time_interpolation
  PUBLIC :: get_current_bc_sst_sic_year

  INTEGER(i8), SAVE :: current_year = -1

  INTEGER                            :: nyears
  INTEGER                            :: imonth_beg, imonth_end

  LOGICAL                            :: lend_of_year

  TYPE(t_time_interpolation_weights) :: tiw_beg
  TYPE(t_time_interpolation_weights) :: tiw_end

CONTAINS
  
  SUBROUTINE read_bc_sst_sic(year, p_patch)
    INTEGER(i8),   INTENT(IN) :: year
    TYPE(t_patch), INTENT(IN) :: p_patch
    INTEGER :: jg
    CHARACTER(len=16) :: fn

    jg = p_patch%id
    lend_of_year = ( time_config%tc_stopdate%date%month  == 1  .AND. &
      &              time_config%tc_stopdate%date%day    == 1  .AND. &
      &              time_config%tc_stopdate%time%hour   == 0  .AND. &
      &              time_config%tc_stopdate%time%minute == 0  .AND. &
      &              time_config%tc_stopdate%time%second == 0 ) 

    nyears = time_config%tc_stopdate%date%year - time_config%tc_startdate%date%year + 1
    IF ( lend_of_year ) nyears = nyears - 1

    ! ----------------------------------------------------------------------

    tiw_beg = calculate_time_interpolation_weights(time_config%tc_startdate)
    tiw_end = calculate_time_interpolation_weights(time_config%tc_stopdate)

    IF ( nyears > 1 ) THEN
      imonth_beg = 0
      imonth_end = 13  
      IF ( tiw_beg%month1_index == 0 ) imonth_end = tiw_end%month2_index
    ELSE
      imonth_beg = tiw_beg%month1_index 
      imonth_end = tiw_end%month2_index
    ENDIF
 
    IF ( lend_of_year ) imonth_end = 13

    WRITE(message_text,'(a,i2,a,i2)') &
       & ' Allocating SST and SIC for months ', imonth_beg, ' to ', imonth_end
    CALL message('mo_bc_sst_sic:read_bc_sst_sic', message_text)

    IF ( imonth_beg > imonth_end ) THEN
      WRITE (message_text, '(a)') 'imonth_beg < imonth_end' 
      CALL finish('mo_bc_sst_sic:read_bc_sst_sic', message_text)
    ENDIF

    IF (n_dom > 1) THEN
      WRITE(fn, '(a,i2.2,a)') 'bc_sst_DOM', jg, '.nc'
    ELSE
      fn = 'bc_sst.nc'
    END IF
    IF (my_process_is_mpi_workroot()) THEN
      WRITE(message_text,'(3a,i0)') 'Read SST from ', TRIM(fn), ' for ', year
      CALL message('',message_text)
    ENDIF
    IF (.NOT.ASSOCIATED(ext_sea(jg)%sst)) THEN
      ALLOCATE (ext_sea(jg)%sst(nproma, p_patch%nblks_c, imonth_beg:imonth_end))
      !$ACC ENTER DATA PCREATE(ext_sea(jg)%sst)
    ENDIF
    CALL read_sst_sic_data(p_patch, ext_sea(jg)%sst, TRIM(fn), year)

    IF (n_dom > 1) THEN
      WRITE(fn, '(a,i2.2,a)') 'bc_sic_DOM', jg, '.nc'
    ELSE
      fn = 'bc_sic.nc'
    END IF
    IF (my_process_is_mpi_workroot()) THEN
      WRITE(message_text,'(3a,i0)') 'Read sea ice from ',TRIM(fn),' for ',year
      CALL message('',message_text)
    ENDIF
    IF (.NOT.ASSOCIATED(ext_sea(jg)%sic)) THEN
      ALLOCATE (ext_sea(jg)%sic(nproma, p_patch%nblks_c, imonth_beg:imonth_end))
      !$ACC ENTER DATA PCREATE(ext_sea(jg)%sic)
    ENDIF
    CALL read_sst_sic_data(p_patch, ext_sea(jg)%sic, TRIM(fn), year)
    
    IF (jg==n_dom) current_year = year

    !$ACC UPDATE DEVICE(ext_sea(jg)%sst, ext_sea(jg)%sic) ASYNC(1)

  END SUBROUTINE read_bc_sst_sic
  
  SUBROUTINE read_sst_sic_data(p_patch, dst, fn, y)
!TODO: switch to reading via mo_read_netcdf_distributed?
    TYPE(t_patch), INTENT(in) :: p_patch
    REAL(dp), CONTIGUOUS, INTENT(INOUT) :: dst(:,:,imonth_beg:)
    CHARACTER(len=*), INTENT(IN) :: fn
    INTEGER(i8), INTENT(in) :: y
    REAL(dp), ALLOCATABLE :: zin(:)
    REAL(dp) :: dummy(0)
    INTEGER :: vlID, taxID, tsID, ts_idx, strID, nmiss, vd, vy, vm, ts_found
    LOGICAL :: found_last_ts, lexist
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: cdiErrorText
    CHARACTER(len=*), PARAMETER :: routine = 'mo_bc_sst_sic:read_sst_sic_data'

    IF (my_process_is_mpi_workroot()) THEN
      INQUIRE (file=fn, exist=lexist)
      IF (.NOT.lexist) THEN
        WRITE (message_text, '(3a)') 'Could not open file ', fn, ': run terminated.'
        CALL finish(routine, message_text)
      ENDIF
      strID = streamOpenRead(fn)
      IF ( strID < 0 ) THEN
        CALL cdiGetStringError(strID, cdiErrorText)
        WRITE (message_text, '(4a)') 'Could not open file ', fn, ': ', cdiErrorText
        CALL finish(routine, message_text)
      END IF
      vlID = streamInqVlist(strID)
      taxID = vlistInqTaxis(vlID)
      tsID = 0
      found_last_ts = .FALSE.
      ts_found = 0

      ALLOCATE(zin(p_patch%n_patch_cells_g))
      DO WHILE (.NOT. found_last_ts)
        IF (streamInqTimestep(strID, tsID) == 0) EXIT
        vd = taxisInqVdate(taxID)
        vy = vd/10000
        vm = (vd/100)-vy*100
        ts_idx = -1

        IF (INT(vy,i8) == y-1_i8 .AND. vm == 12) THEN
          IF ( imonth_beg == 0 ) THEN
              ts_idx = 0
              ts_found = ts_found + 1
          END IF
        ELSE IF (INT(vy,i8) == y) THEN
          IF ( vm >= imonth_beg .AND. vm <= imonth_end ) THEN 
            ts_idx = vm
            ts_found = ts_found + 1
            IF ( vm == imonth_end ) found_last_ts = .TRUE.
          END IF
        ELSE IF (INT(vy,i8) == y+1_i8 .AND. vm == 1) THEN
          IF ( imonth_end == 13 ) THEN
            ts_idx = 13
            ts_found = ts_found + 1
            found_last_ts = .TRUE.
          END IF
        END IF
        IF (ts_idx /= -1) THEN
          CALL streamReadVarSlice(strID, 0, 0, zin, nmiss)
          CALL p_bcast(ts_idx, process_mpi_root_id, p_comm_work)
          dst(:,SIZE(dst,2),ts_idx) = 0._dp
          CALL p_patch%comm_pat_scatter_c%distribute(zin, dst(:,:,ts_idx), .FALSE.)
        ENDIF
        tsID = tsID+1
      END DO
      IF (ts_found < imonth_end - imonth_beg + 1) &
          & CALL finish ('mo_bc_sst_sic:read_sst_sic_data', &
          &   'could not read required data from input file')
      ts_idx = -1
      CALL p_bcast(ts_idx, process_mpi_root_id, p_comm_work)
      DEALLOCATE(zin)
      CALL streamClose(strID)
    ELSE
      ts_idx = 0
      DO
        CALL p_bcast(ts_idx, process_mpi_root_id, p_comm_work)
        IF(ts_idx .EQ. -1) EXIT
        dst(:,SIZE(dst,2),ts_idx) = 0._dp
        CALL p_patch%comm_pat_scatter_c%distribute(dummy, dst(:,:,ts_idx), .FALSE.)
      END DO
    END IF
  END SUBROUTINE read_sst_sic_data

  SUBROUTINE bc_sst_sic_time_interpolation(tiw, tsw, seaice, siced, p_patch, mask, l_init, lopenacc)
    
    TYPE( t_time_interpolation_weights), INTENT(in) :: tiw
    REAL(dp)       , INTENT(inout) :: tsw(:,:) 
    REAL(dp)       , INTENT(out) :: seaice(:,:) 
    REAL(dp)       , INTENT(out) :: siced(:,:) 
    TYPE(t_patch)  , INTENT(in)  :: p_patch
    LOGICAL        , INTENT(in)  :: mask(:,:)  !< logical mask, indicating where to apply tsw and sea ice/depth
    LOGICAL        , INTENT(in)  :: l_init     !< switch for first call at initialization
    LOGICAL, OPTIONAL, INTENT(in) :: lopenacc  ! Flag to run on GPU
    ! If l_init=.FALSE., tsw is only computed where mask==.TRUE (this is used
    ! at all time steps in the time loop). If l_init=.TRUE., tsw is computed 
    ! everywhere (this is used during initialization in order to initialize ts_tile(:,:,iwtr)
    ! also over land/lakes). 
    ! Note that lakes and ocean/sea ice are mutually exclusive, i.e. a cell cannot
    ! contain both lake and ocean/sea ice.

    REAL(dp) :: zts(SIZE(tsw,1),SIZE(tsw,2))
    REAL(dp) :: zic(SIZE(tsw,1),SIZE(tsw,2))
    REAL(dp) :: ztsw(SIZE(tsw,1),SIZE(tsw,2))
    REAL(dp), CONTIGUOUS, POINTER :: sst(:,:,:), sic(:,:,:)

    INTEGER  :: jc, jb, jg, jce, nblk

    jg = p_patch%id

    sst => ext_sea(jg)%sst
    sic => ext_sea(jg)%sic

    jce  = SIZE(tsw,1)
    nblk = SIZE(tsw,2)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) COPYIN(tiw) CREATE(zts, zic, ztsw) IF(lzopenacc)
!$omp parallel private(jb,jc)
!$omp do
    !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
    DO jb = 1, nblk
      DO jc = 1, jce
        zts(jc,jb) = tiw%weight1 * sst(jc,jb,tiw%month1_index) + tiw%weight2 * sst(jc,jb,tiw%month2_index)
        zic(jc,jb) = tiw%weight1 * sic(jc,jb,tiw%month1_index) + tiw%weight2 * sic(jc,jb,tiw%month2_index)
      END DO
    END DO
!$omp end do nowait

    !TODO: missing siced needs to be added

!$omp do
    !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
    DO jb = 1, nblk
      DO jc = 1, jce
        ! assuming input data is in percent
        seaice(jc,jb) = zic(jc,jb) * MERGE(0.01_dp, 0.0_dp, mask(jc,jb))
      END DO
    END DO
!$omp end do nowait

!$omp do
    !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
    DO jb = 1, nblk
      DO jc = 1, jce
        seaice(jc,jb) = MERGE(0.99_dp, seaice(jc,jb), seaice(jc,jb) > 0.99_dp)
        seaice(jc,jb) = MERGE(0.0_dp, seaice(jc,jb), seaice(jc,jb) <= 0.01_dp)

        ztsw(jc,jb) = MERGE(tf_salt, MAX(zts(jc,jb), tf_salt), seaice(jc,jb) > 0.0_dp) 
      END DO
    END DO
!$omp end do nowait

    IF (l_init) THEN
!$omp do
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO jb = 1, nblk
        DO jc = 1, jce
          tsw(jc,jb) = ztsw(jc,jb)
        END DO
      END DO
!$omp end do nowait
    ELSE
!$omp do
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO jb = 1, nblk
        DO jc = 1, jce
          IF (mask(jc,jb)) THEN
            tsw(jc,jb) = ztsw(jc,jb)
          END IF
        END DO
      END DO
!$omp end do nowait
    END IF

!$omp do
    !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
    DO jb = 1, nblk
      DO jc = 1, jce
        IF (seaice(jc,jb) > 0.0_dp) THEN
          siced(jc,jb) = MERGE(2._dp, 1._dp, p_patch%cells%center(jc,jb)%lat > 0.0_dp)
        ELSE
          siced(jc,jb) = 0._dp
        END IF
      END DO
    END DO
!$omp end do nowait
!$omp end parallel
    !$ACC END PARALLEL
    !$ACC WAIT(1)

    !CALL message('','Interpolated sea surface temperature and sea ice cover.')

  END SUBROUTINE bc_sst_sic_time_interpolation

  FUNCTION get_current_bc_sst_sic_year() RESULT(this_year)
    INTEGER(i8) :: this_year
    this_year = current_year
  END FUNCTION get_current_bc_sst_sic_year

END MODULE mo_bc_sst_sic
