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

! Subroutine read_bc_ozone reads monthly ozone
! concentrations from yearly files. The ozone concentrations are
! used as boundary conditions for the radiative forcing of the
! atmosphereAMIP. The routine is called from mo_aes_phy_interface.f90.

MODULE mo_bc_ozone

  USE mo_kind,                     ONLY: wp, i8
  USE mo_exception,                ONLY: message, message_text, warning
  USE mo_model_domain,             ONLY: t_patch
  USE mo_parallel_config,          ONLY: p_test_run, nproma
  USE mo_io_config,                ONLY: default_read_method
  USE mo_read_interface,           ONLY: openInputFile, closeFile,        &
  &                                      read_3D_time, t_stream_id, on_cells
  USE mo_mpi,                      ONLY: my_process_is_stdio, p_bcast,        &
  &                                      p_comm_work_test, p_comm_work, p_io
  USE mo_impl_constants,           ONLY: max_dom
  USE mo_grid_config,              ONLY: n_dom
  USE mo_aes_rad_config,           ONLY: aes_rad_config
  USE mo_netcdf_errhandler,        ONLY: nf
  USE mo_netcdf
  USE mo_time_config,              ONLY: time_config
  USE mo_bcs_time_interpolation,   ONLY: t_time_interpolation_weights, &
       &                                      calculate_time_interpolation_weights

  USE mo_timer,                     ONLY: ltimer, timer_start, timer_stop, &
                                          timer_coupling
  USE mo_atmo_o3_provider_coupling, ONLY: nplev_o3_provider, &
                                          couple_atmo_to_o3_provider

  IMPLICIT NONE
  PRIVATE
  INTEGER(i8), SAVE                 :: pre_year(max_dom)=-HUGE(1) ! Variable to check if it is time to read

  PUBLIC                            :: ext_ozone
  PUBLIC                            :: read_bc_ozone

  TYPE t_ext_ozone
    REAL(wp), ALLOCATABLE           :: o3_plev(:,:,:,:)     ! Monthly ozone mass mixing ratio at pressure levels
    REAL(wp), ALLOCATABLE           :: plev_full_o3(:)      ! Full pressure levels
    REAL(wp), ALLOCATABLE           :: plev_half_o3(:)      ! Half pressure levels, derived from plev_full_o3
    INTEGER                         :: nplev_o3             ! Number of full pressure levels
    INTEGER(i8)                     :: year                 ! Reference year.
  END TYPE t_ext_ozone

  TYPE(t_ext_ozone), ALLOCATABLE, TARGET :: ext_ozone(:)

  TYPE(t_time_interpolation_weights) :: tiw_beg, tiw_end
  INTEGER                            :: nyears
  LOGICAL                            :: lend_of_year

CONTAINS

  SUBROUTINE read_bc_ozone(year, p_patch, irad_o3, vmr2mmr_opt, opt_from_coupler)

    INTEGER(i8)  , INTENT(in)            :: year
    TYPE(t_patch), TARGET, INTENT(in)    :: p_patch
    INTEGER      , INTENT(in)            :: irad_o3
    REAL(wp)     , INTENT(in), OPTIONAL  :: vmr2mmr_opt
    LOGICAL      , INTENT(IN), OPTIONAL  :: opt_from_coupler

    CHARACTER(len=512)                :: fname
    TYPE(t_stream_id)                 :: stream_id
    CHARACTER(len=4)                  :: cyear
    CHARACTER(len=2)                  :: cjg
    CHARACTER(len=*), PARAMETER       :: subprog_name = 'mo_bc_ozone:read_bc_ozone'
    INTEGER                           :: ncid, varid, mpi_comm, jg
    INTEGER                           :: nplev_o3
    INTEGER                           :: imonth_beg, imonth_end      ! months with range 0-13
    INTEGER                           :: kmonth_beg, kmonth_end      ! months with range 1-12
    INTEGER                           :: nmonths
    REAL(wp), POINTER                 :: zo3_plev(:,:,:,:)           ! (nproma, levels, blocks, time)
    REAL(wp)                          :: vmr2mmr_o3
    LOGICAL                           :: l_first
    LOGICAL                           :: from_coupler = .FALSE.

    jg    = p_patch%id
    WRITE(cjg,'(i2.2)') jg

    IF ( PRESENT(vmr2mmr_opt) ) THEN
      vmr2mmr_o3 = vmr2mmr_opt
    ELSE
      vmr2mmr_o3 = 1._wp
    END IF

    ! This should be part of an initialisation routine
    ! check when mergeing with branch icon-aes-link-echam-bc

    tiw_beg = calculate_time_interpolation_weights(time_config%tc_startdate)
    tiw_end = calculate_time_interpolation_weights(time_config%tc_stopdate)

    lend_of_year = ( time_config%tc_stopdate%date%month  == 1  .AND. &
      &              time_config%tc_stopdate%date%day    == 1  .AND. &
      &              time_config%tc_stopdate%time%hour   == 0  .AND. &
      &              time_config%tc_stopdate%time%minute == 0  .AND. &
      &              time_config%tc_stopdate%time%second == 0 ) 

    nyears = time_config%tc_stopdate%date%year - time_config%tc_startdate%date%year + 1
    IF ( lend_of_year ) nyears = nyears - 1

    ! ------------------------------------------------

    ! allocate once only structure for all grids
    IF (.NOT. ALLOCATED(ext_ozone)) ALLOCATE(ext_ozone(n_dom))
    !$ACC ENTER DATA PCREATE(ext_ozone)

    l_first = .NOT. ALLOCATED(ext_ozone(jg)% o3_plev)

    IF (PRESENT(opt_from_coupler)) from_coupler = opt_from_coupler

    IF (from_coupler) THEN

       nplev_o3 =  nplev_o3_provider
       ext_ozone(jg)%nplev_o3 = nplev_o3_provider
       IF ( .NOT. ALLOCATED(ext_ozone(jg)% o3_plev) ) THEN
          ALLOCATE(ext_ozone(jg)% o3_plev(nproma, nplev_o3, p_patch%nblks_c, 1))
          ext_ozone(jg)% o3_plev = 0.0_wp
       END IF

       IF (ltimer) CALL timer_start(timer_coupling)
       CALL couple_atmo_to_o3_provider( &
         p_patch, vmr2mmr_o3, ext_ozone(jg)% o3_plev)
       IF (ltimer) CALL timer_stop(timer_coupling)

       fname = 'bc_ozone.nc'

    END IF

    IF (year > pre_year(jg) .AND. .NOT. from_coupler) THEN
      !
      ! If year = pre_year, then the external monthly ozone data are already stored.
      ! Nothing needs to be done.
      !
      IF (ALLOCATED(ext_ozone(jg)% o3_plev)) THEN
        !
        ! If the memory for o3_plev is already allocated, then this is a new
        ! year of a simulation that started before this year, meaning that
        ! a year of external monthly ozone data is already stored.
        !

        IF (irad_o3==5) THEN
          !
          ! ozone is transient and the external monthly ozone data
          ! of this year and January of the next year must be read from file.
          !
          ! For other irad_o3 cases no new data must be read.
          ! Nothing needs to be done.
          !
          WRITE(message_text,'(2a)') 'Copy ozone for months 12:13 to months 0:1'
          CALL message('read_bc_ozone', message_text)
          ext_ozone(jg)% o3_plev(:,:,:,0:1) = ext_ozone(jg)% o3_plev(:,:,:,12:13)

          WRITE(cyear,'(i4)') year
          IF (n_dom > 1) THEN
            fname = 'bc_ozone_DOM'//TRIM(cjg)//'_'//TRIM(cyear)//'.nc'
          ELSE
            fname = 'bc_ozone_'//TRIM(cyear)//'.nc'
          END IF

          IF ( year < time_config%tc_stopdate%date%year ) THEN
            imonth_end = 13
          ELSE
            imonth_end = tiw_end%month2_index
          ENDIF

          kmonth_end = MIN(imonth_end,12)

          IF ( kmonth_end > 1 ) THEN
            WRITE(message_text,'(2a)') 'Read ozone for months  2:12 from file: ',TRIM(fname)
            CALL message('read_bc_ozone', message_text)
            CALL openInputFile(stream_id, fname, p_patch, default_read_method)
            CALL read_3D_time(stream_id=stream_id, location=on_cells,      &
              &               variable_name='O3', return_pointer=zo3_plev, &
              &               start_timestep=2,end_timestep=kmonth_end)
            CALL closeFile(stream_id)
            ext_ozone(jg)% o3_plev(:,:,:,2:kmonth_end) = vmr2mmr_o3*zo3_plev(:,:,:,1:kmonth_end-1)

            WRITE(cyear,'(i4)') year+1
            IF (n_dom > 1) THEN
              fname = 'bc_ozone_DOM'//TRIM(cjg)//'_'//TRIM(cyear)//'.nc'
            ELSE
              fname = 'bc_ozone_'//TRIM(cyear)//'.nc'
            END IF
          ENDIF

          IF ( imonth_end == 13 ) THEN
            WRITE(message_text,'(2a)') 'Read ozone for month     13 from file: ',TRIM(fname)
            CALL message('read_bc_ozone', message_text)
            CALL openInputFile(stream_id, fname, p_patch, default_read_method)
            CALL read_3D_time(stream_id=stream_id, location=on_cells,      &
              &               variable_name='O3', return_pointer=zo3_plev, &
              &               start_timestep=1,end_timestep=1)
            CALL closeFile(stream_id)
            ext_ozone(jg)% o3_plev(:,:,:,imonth_end) = vmr2mmr_o3*zo3_plev(:,:,:,1)
          ENDIF

          ext_ozone(jg)%year = year
        END IF
        !
      ELSE
        !
        ! If the memory of o3_plev is not allocated, then this is
        ! a new simulation or a simulation that starts from a restart file
        ! and ozone data must be read.
        !
        ! Depending on the irad_o3 case different amounts of data must be read.
        !
        SELECT CASE (irad_o3)
        !
        CASE (6) ! Ozone has a climatological annual cycle defined by monthly data in an annual file
          !
          IF ( nyears > 1 ) THEN
            imonth_beg = 0
            imonth_end = 13
          ELSE
            imonth_beg = tiw_beg%month1_index
            imonth_end = tiw_end%month2_index
            IF ( lend_of_year ) imonth_end = 13
          ENDIF
          !
          kmonth_beg=MAX( 1,imonth_beg)
          kmonth_end=MIN(12,imonth_end)
          nmonths = kmonth_end-kmonth_beg+1
          !
          IF (n_dom > 1) THEN
            fname = 'bc_ozone_DOM'//TRIM(cjg)//'.nc'
          ELSE
            fname = 'bc_ozone'//'.nc'
          END IF
           !
          CALL openInputFile(stream_id, fname, p_patch, default_read_method)
          !
          WRITE(message_text,'(a,i2.2,a1,i2.2,2a)') 'Read clim. annual cycle of ozone for months ', &
            &                                        kmonth_beg, ':', kmonth_end, ' from file: ', fname
          CALL message('read_bc_ozone', message_text)
          !
          CALL read_3D_time(stream_id=stream_id, location=on_cells,      &
            &               variable_name='O3', return_pointer=zo3_plev, &
            &               start_timestep=kmonth_beg,end_timestep=kmonth_end)
          !
          ! Now the spatial dimensions are known --> allocate memory for months 0:13
          WRITE(message_text,'(a,i2.2,a1,i2.2)') 'Alloc clim. annual cycle of ozone for months ', &
            &                                        imonth_beg, ':', imonth_end
          CALL message('read_bc_ozone', message_text)
          ALLOCATE(ext_ozone(jg)% o3_plev(SIZE(zo3_plev,1), &
            &                             SIZE(zo3_plev,2), &
            &                             SIZE(zo3_plev,3),0:13))
          !$ACC ENTER DATA PCREATE(ext_ozone(jg)%o3_plev)
          ext_ozone(jg)% o3_plev = 0.0_wp
          !
          ext_ozone(jg)% o3_plev(:,:,:,kmonth_beg:kmonth_end) = vmr2mmr_o3*zo3_plev(:,:,:,1:nmonths)
          !
          IF ( imonth_beg == 0 ) THEN
            IF ( kmonth_end < 12 ) THEN
              !
              WRITE(message_text,'(2a)') 'Read clim. annual cycle of ozone for month 12 ( -> 0) from file: ', TRIM(fname)
              CALL message('read_bc_ozone', message_text)
              !
              CALL read_3D_time(stream_id=stream_id, location=on_cells,      &
                &               variable_name='O3', return_pointer=zo3_plev, &
                &               start_timestep=12,end_timestep=12)
              ext_ozone(jg)% o3_plev(:,:,:,0) = vmr2mmr_o3*zo3_plev(:,:,:,1)
            ELSE
              ext_ozone(jg)% o3_plev(:,:,:,0) = vmr2mmr_o3*zo3_plev(:,:,:,12)
            ENDIF
          ENDIF
          !
          IF ( imonth_end == 13 ) THEN
            IF ( kmonth_beg > 1 ) THEN
              !
              WRITE(message_text,'(2a)') 'Read clim. annual cycle of ozone for month 1 ( -> 13) from file: ', TRIM(fname)
              CALL message('read_bc_ozone', message_text)
              !
              CALL read_3D_time(stream_id=stream_id, location=on_cells,      &
                &               variable_name='O3', return_pointer=zo3_plev, &
                &               start_timestep=1,end_timestep=1)
            ENDIF
            ext_ozone(jg)% o3_plev(:,:,:,13) = vmr2mmr_o3*zo3_plev(:,:,:,1)
          ENDIF
          !
          CALL closeFile(stream_id)
          !
        CASE (4) ! Ozone is constant in time
          !
          IF (n_dom > 1) THEN
            fname = 'bc_ozone_DOM'//TRIM(cjg)//'.nc'
          ELSE
            fname = 'bc_ozone'//'.nc'
          END IF
          !
          WRITE(message_text,'(2a)') 'Read constant-in-time ozone from file: ',TRIM(fname)
          CALL message('read_bc_ozone', message_text)
          !
          CALL openInputFile(stream_id, fname, p_patch, default_read_method)
          CALL read_3D_time(stream_id=stream_id, location=on_cells,      &
            &               variable_name='O3', return_pointer=zo3_plev, &
            &               start_timestep=1,end_timestep=1)
          CALL closeFile(stream_id)
          !
          ! Now the spatial dimensions are known --> allocate memory for one time slice
          ALLOCATE(ext_ozone(jg)% o3_plev(SIZE(zo3_plev,1),SIZE(zo3_plev,2),SIZE(zo3_plev,3),1))
          ext_ozone(jg)% o3_plev = 0.0_wp
          !$ACC ENTER DATA PCREATE(ext_ozone(jg)%o3_plev)
          !
          ext_ozone(jg)% o3_plev(:,:,:,1) = vmr2mmr_o3*zo3_plev(:,:,:,1)
          !
          !

        CASE (5) ! Ozone is transient and defined by monthly data in annual files

          !
          IF ( nyears > 1 ) THEN
            imonth_beg = 0
            imonth_end = 13
          ELSE
            imonth_beg = tiw_beg%month1_index
            imonth_end = tiw_end%month2_index
            IF ( lend_of_year ) imonth_end = 13
          ENDIF
          !
          kmonth_beg=MAX( 1,imonth_beg)
          kmonth_end=MIN(12,imonth_end)
          nmonths = kmonth_end-kmonth_beg+1
          !
          ! 1. Read December of the previous year
          !
          IF ( imonth_beg == 0 ) THEN
            !
            WRITE(cyear,'(i4)') year-1
            !
            IF (n_dom > 1) THEN
              fname = 'bc_ozone_DOM'//TRIM(cjg)//'_'//TRIM(cyear)//'.nc'
            ELSE
              fname = 'bc_ozone_'//TRIM(cyear)//'.nc'
            END IF
            !
            WRITE(message_text,'(2a)') 'Read ozone for month      0 from file: ',TRIM(fname)
            CALL message('read_bc_ozone', message_text)
            CALL openInputFile(stream_id, fname, p_patch, default_read_method)
            CALL read_3D_time(stream_id=stream_id, location=on_cells,         &
              &               variable_name='O3', return_pointer=zo3_plev, &
              &               start_timestep=12,end_timestep=12)
            CALL closeFile(stream_id)
            !
            ! Now the spatial dimensions are known -->
            !       allocate memory for months imonth_beg to imonth_end
            ! o3_plev has to be allocated from 0:13, because
            ! time dimensions are hardcoded in o3_timeint as well
            ALLOCATE(ext_ozone(jg)% o3_plev(SIZE(zo3_plev,1), &
                                            SIZE(zo3_plev,2), &
                                            SIZE(zo3_plev,3), 0:13))
            ext_ozone(jg)% o3_plev = 0.0_wp
            !$ACC ENTER DATA PCREATE(ext_ozone(jg)%o3_plev)
            !
            ext_ozone(jg)% o3_plev(:,:,:,0) = vmr2mmr_o3*zo3_plev(:,:,:,1)
          ENDIF
          !
          ! 2. Read January-December of this year
          !
          WRITE(cyear,'(i4)') year
          IF (n_dom > 1) THEN
            fname = 'bc_ozone_DOM'//TRIM(cjg)//'_'//TRIM(cyear)//'.nc'
          ELSE
            fname = 'bc_ozone_'//TRIM(cyear)//'.nc'
          END IF
          !
          WRITE(message_text,'(a,i2.2,a1,i2.2,2a)') 'Read ozone for months ', kmonth_beg, ':', &
          &                            kmonth_end, ' from file: ',TRIM(fname)
          CALL message('read_bc_ozone', message_text)
          CALL openInputFile(stream_id, fname, p_patch, default_read_method)
          CALL read_3D_time(stream_id=stream_id, location=on_cells,         &
            &               variable_name='O3', return_pointer=zo3_plev, &
            &               start_timestep=kmonth_beg,end_timestep=kmonth_end)
          CALL closeFile(stream_id)
          !
          ! Now the spatial dimensions are known -->
          !       allocate memory for months imonth_beg to imonth_end
          IF ( imonth_beg > 0 ) THEN
            ! o3_plev has to be allocated from 0:13, because
            ! time dimensions are hardcoded in o3_timeint as well
            ALLOCATE(ext_ozone(jg)% o3_plev(SIZE(zo3_plev,1), &
                                            SIZE(zo3_plev,2), &
                                            SIZE(zo3_plev,3), 0:13))
            ext_ozone(jg)% o3_plev = 0.0_wp
            !$ACC ENTER DATA PCREATE(ext_ozone(jg)%o3_plev)
          ENDIF
          !
          ext_ozone(jg)% o3_plev(:,:,:,kmonth_beg:kmonth_end) = vmr2mmr_o3*zo3_plev(:,:,:,1:nmonths)
          !
          ! 3. Read January of the next year
          !
          IF ( imonth_end == 13 ) THEN
            WRITE(cyear,'(i4)') year+1
            IF (n_dom > 1) THEN
              fname = 'bc_ozone_DOM'//TRIM(cjg)//'_'//TRIM(cyear)//'.nc'
            ELSE
              fname = 'bc_ozone_'//TRIM(cyear)//'.nc'
            END IF
            !
            WRITE(message_text,'(2a)') 'Read ozone for month     13 from file: ',TRIM(fname)
            CALL message('read_bc_ozone', message_text)
            CALL openInputFile(stream_id, fname, p_patch, default_read_method)
            CALL read_3D_time(stream_id=stream_id, location=on_cells,         &
               &              variable_name='O3', return_pointer=zo3_plev, &
               &              start_timestep=1,end_timestep=1)
            CALL closeFile(stream_id)
            !
            ! Safety dance, as this should not happen, but who knows ...
            ! Now the spatial dimensions are known -->
            !       allocate memory for month 13 only
            IF ( imonth_beg == 13 ) THEN
              ! o3_plev has to be allocated from 0:13, because
              ! time dimensions are hardcoded in o3_timeint as well
              ALLOCATE(ext_ozone(jg)% o3_plev(SIZE(zo3_plev,1), &
                                              SIZE(zo3_plev,2), &
                                              SIZE(zo3_plev,3), 0:13))
              ext_ozone(jg)% o3_plev = 0.0_wp
              !$ACC ENTER DATA PCREATE(ext_ozone(jg)%o3_plev)
            ENDIF
            !
            ext_ozone(jg)% o3_plev(:,:,:,13) = vmr2mmr_o3*zo3_plev(:,:,:,1)
            !
          ENDIF

          ext_ozone(jg)%year = year

        END SELECT

      END IF
      !$ACC UPDATE DEVICE(ext_ozone(jg)%o3_plev) ASYNC(1)

      pre_year(jg) = year

    END IF ! (year > pre_year(jg))  .AND. .NOT. from_coupler

    IF (l_first) THEN

      ! Read pressure level grid on which the external ozone data are valid
      ! This needs to be done equally for all irad_o3 cases.

      nplev_o3 = SIZE(ext_ozone(jg)% o3_plev,2)

      ext_ozone(jg)% nplev_o3 = nplev_o3

      !$ACC WAIT(1)
      !$ACC EXIT DATA DELETE(ext_ozone(jg)%plev_full_o3) IF(ALLOCATED(ext_ozone(jg)%plev_full_o3))
      !$ACC EXIT DATA DELETE(ext_ozone(jg)%plev_half_o3) IF(ALLOCATED(ext_ozone(jg)%plev_half_o3))
      IF(ALLOCATED(ext_ozone(jg)% plev_full_o3)) DEALLOCATE(ext_ozone(jg)% plev_full_o3)
      IF(ALLOCATED(ext_ozone(jg)% plev_half_o3)) DEALLOCATE(ext_ozone(jg)% plev_half_o3)
      ALLOCATE(ext_ozone(jg)% plev_full_o3(nplev_o3  ))
      ALLOCATE(ext_ozone(jg)% plev_half_o3(nplev_o3+1))
      !$ACC ENTER DATA PCREATE(ext_ozone(jg)%plev_full_o3, ext_ozone(jg)%plev_half_o3)

      mpi_comm = MERGE(p_comm_work_test, p_comm_work, p_test_run)

      IF(my_process_is_stdio()) THEN
        CALL nf(nf90_open(TRIM(fname), NF90_NOWRITE, ncid), subprog_name)
        CALL nf(nf90_inq_varid(ncid, 'plev', varid), subprog_name)
        CALL nf(nf90_get_var(ncid, varid, ext_ozone(jg)% plev_full_o3), subprog_name)
        CALL nf(nf90_close(ncid), subprog_name)
      END IF
      CALL p_bcast(ext_ozone(jg)% plev_full_o3(:), p_io, mpi_comm)

      ! define half levels of ozone pressure grid
      ! upper boundary: ph =      0.Pa -> extrapolation of uppermost value
      ! lower boundary: ph = 125000.Pa -> extrapolation of lowermost value
      ext_ozone(jg)% plev_half_o3(1)          = 0._wp
      ext_ozone(jg)% plev_half_o3(2:nplev_o3) = 0.5_wp*( ext_ozone(jg)% plev_full_o3(1:nplev_o3-1) &
        &                                               +ext_ozone(jg)% plev_full_o3(2:nplev_o3  ) )

      ext_ozone(jg)% plev_half_o3(nplev_o3+1) = 125000._wp

      l_first = .FALSE.

      ! Set pointer for OpenACC
      !$ACC UPDATE DEVICE(ext_ozone(jg)%plev_half_o3, ext_ozone(jg)%plev_full_o3) ASYNC(1)

    END IF

  END SUBROUTINE read_bc_ozone

END MODULE mo_bc_ozone
