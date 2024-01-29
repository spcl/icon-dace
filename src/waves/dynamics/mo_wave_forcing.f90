! Reading of forcing data in standalone mode
!
! Routines for reading forcing data at user-defined time intervals in
! standalone mode.
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

MODULE mo_wave_forcing

  USE mo_kind,                     ONLY: wp
  USE mo_exception,                ONLY: finish, warning, message, message_text
  USE mo_io_units,                 ONLY: filename_max
  USE mo_model_domain,             ONLY: t_patch
  USE mo_reader_sst_sic,           ONLY: t_sst_sic_reader
  USE mo_interpolate_time,         ONLY: t_time_intp
  USE mo_fortran_tools,            ONLY: copy, DO_DEALLOCATE
  USE mtime,                       ONLY: datetime
  USE mo_wave_td_update,           ONLY: update_ice_free_mask, &
    &                                    update_speed_and_direction
  USE mo_mpi,                      ONLY: my_process_is_mpi_workroot, p_io, p_bcast, &
    &                                    p_comm_work

  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_forcing'

  ! wave-specific forcing data type
  !
  TYPE :: t_read_wave_forcing

    TYPE(t_patch), POINTER :: p_patch => NULL()  ! pointer to actual patch
    ! time interpolated fields
    !
    REAL(wp), ALLOCATABLE :: u10m_raw (:,:,:,:)  ! zonal wind in 10m (m/s)
    REAL(wp), ALLOCATABLE :: v10m_raw (:,:,:,:)  ! meridional wind in 10m (m/s)
    REAL(wp), ALLOCATABLE :: sic_raw  (:,:,:,:)  ! sea ice concentration at centers (1)
    REAL(wp), ALLOCATABLE :: slh_raw  (:,:,:,:)  ! sea level height at centers ( m )
    REAL(wp), ALLOCATABLE :: uosc_raw (:,:,:,:)  ! zonal ocean surface current at centers ( m/s )
    REAL(wp), ALLOCATABLE :: vosc_raw (:,:,:,:)  ! meridional ocean surface current at centers ( m/s )

    ! reader
    TYPE(t_sst_sic_reader) :: uv10_reader
    TYPE(t_sst_sic_reader) :: sic_reader
    TYPE(t_sst_sic_reader) :: slh_reader
    TYPE(t_sst_sic_reader) :: osc_reader

    ! time interpolator
    TYPE(t_time_intp) :: u10_intp, v10_intp
    TYPE(t_time_intp) :: sic_intp
    TYPE(t_time_intp) :: slh_intp
    TYPE(t_time_intp) :: uosc_intp, vosc_intp

    ! forcing flags
    LOGICAL :: l_wind_exist, l_ice_exist, l_slh_exist, l_osc_exist

  CONTAINS

   PROCEDURE :: init               => read_wave_forcing__init
   PROCEDURE :: deinit             => read_wave_forcing__deinit
   PROCEDURE :: get_new_rawdata    => read_wave_forcing__get_new_rawdata
   PROCEDURE :: update_forcing     => read_wave_forcing__update_forcing
  END TYPE t_read_wave_forcing


  PUBLIC  :: t_read_wave_forcing


CONTAINS

  !>
  !! Initialize reader of wave forcing data
  !!
  !! Initialize the reader object which reads wave forcing data from file
  !! and performs a linear interpolation in time to a specified target date.
  !!
  SUBROUTINE read_wave_forcing__init (self, p_patch, destination_time, wave_forc_wind_file, &
    &        wave_forc_ice_file, wave_forc_slh_file, wave_forc_osc_file)

    CHARACTER(len=*), PARAMETER :: routine = modname//':read_wave_forcing__init'

    CLASS(t_read_wave_forcing),  TARGET     :: self
    TYPE(t_patch),  TARGET,      INTENT(IN) :: p_patch
    TYPE(datetime), POINTER,     INTENT(IN) :: destination_time
    CHARACTER(LEN=filename_max), INTENT(IN) :: wave_forc_wind_file
    CHARACTER(LEN=filename_max), INTENT(IN) :: wave_forc_ice_file
    CHARACTER(LEN=filename_max), INTENT(IN) :: wave_forc_slh_file
    CHARACTER(LEN=filename_max), INTENT(IN) :: wave_forc_osc_file

  !-----------------------------------------

    ! set pointer to current patch
    self%p_patch => p_patch


    IF(my_process_is_mpi_workroot()) THEN
      INQUIRE (FILE=wave_forc_wind_file, EXIST=self%l_wind_exist)
      IF (self%l_wind_exist) THEN
        CALL message(routine,'10m wind from: '//wave_forc_wind_file)
      ELSE
        WRITE(message_text,'(a,a,a)') 'Instant wind data file ', TRIM(wave_forc_wind_file), &
          &                           ' is not found.'
        CALL finish(routine, message_text)
      ENDIF

      INQUIRE (FILE=wave_forc_ice_file, EXIST=self%l_ice_exist)
      IF (self%l_ice_exist) THEN
        CALL message(routine,'ice concentration from: '//wave_forc_ice_file)
      ELSE
        WRITE(message_text,'(a,a,a)') 'Instant sea-ice data file ', TRIM(wave_forc_ice_file), &
          &                           ' is not found. Run without sea ice'
        CALL warning(routine, message_text)
      ENDIF

      INQUIRE (FILE=wave_forc_slh_file, EXIST=self%l_slh_exist)
      IF (self%l_slh_exist) THEN
        CALL message(routine,'sea level height from: '//wave_forc_slh_file)
      ELSE
        WRITE(message_text,'(a,a,a)') 'Instant sea level data file ', TRIM(wave_forc_slh_file), &
          &                           ' is not found. Run without sea level height'
        CALL warning(routine, message_text)
      ENDIF

      INQUIRE (FILE=wave_forc_osc_file, EXIST=self%l_osc_exist)
      IF (self%l_osc_exist) THEN
        CALL message(routine,'ocean surface currents from: '//wave_forc_osc_file)
      ELSE
        WRITE(message_text,'(a,a,a)') 'Instant ocean surface current data file ', TRIM(wave_forc_osc_file), &
          &                           ' is not found. Run without ocean currents'
        CALL warning(routine, message_text)
      ENDIF
    ENDIF
    !
    ! broadcast l_wind_exist, l_ice_exist, l_slh_exist, l_osc_exist from I-PE to WORK PEs
    CALL p_bcast(self%l_wind_exist, p_io, p_comm_work)
    CALL p_bcast(self%l_ice_exist,  p_io, p_comm_work)
    CALL p_bcast(self%l_slh_exist,  p_io, p_comm_work)
    CALL p_bcast(self%l_osc_exist,  p_io, p_comm_work)



    ! Initialize reader and time interpolator
    IF (self%l_wind_exist) THEN
      CALL self%uv10_reader%init(p_patch, wave_forc_wind_file)
      CALL self%u10_intp   %init(self%uv10_reader, destination_time, "u_10m")
      CALL self%v10_intp   %init(self%uv10_reader, destination_time, "v_10m")
    END IF

    IF (self%l_ice_exist) THEN
      CALL self%sic_reader%init(p_patch, wave_forc_ice_file)
      CALL self%sic_intp  %init(self%sic_reader, destination_time, "fr_seaice")
    END IF

    IF (self%l_slh_exist) THEN
      CALL self%slh_reader%init(p_patch, wave_forc_slh_file)
      CALL self%slh_intp  %init(self%slh_reader, destination_time, "SLH")
    END IF

    IF (self%l_osc_exist) THEN
      CALL self%osc_reader%init(p_patch, wave_forc_osc_file)
      CALL self%uosc_intp  %init(self%osc_reader, destination_time, "UOSC")
      CALL self%vosc_intp  %init(self%osc_reader, destination_time, "VOSC")
    END IF

  END SUBROUTINE read_wave_forcing__init


  !>
  !! Read new timelevel of wave forcing data
  !!
  !! A new timelevel of wave forcing data is read from file and
  !! a linear interpolation to the specified destination time is performed.
  !! The result is stored inside the read object itself.
  !!
  SUBROUTINE read_wave_forcing__get_new_rawdata (self, destination_time)

    CHARACTER(len=*), PARAMETER :: routine = modname//':read_wave_forcing__get_new_rawdata'

    CLASS(t_read_wave_forcing),   TARGET    :: self
    TYPE(datetime), POINTER,     INTENT(IN) :: destination_time  ! validity time of the new data

    IF (self%l_wind_exist) THEN
      ! get new u10m
      CALL self%u10_intp%intp(destination_time, self%u10m_raw)
      ! get new v10m
      CALL self%v10_intp%intp(destination_time, self%v10m_raw)
    END IF

    ! get new sic
    IF (self%l_ice_exist) THEN
      CALL self%sic_intp%intp(destination_time, self%sic_raw)
    ENDIF

    ! get new slh
    IF (self%l_slh_exist) THEN
      CALL self%slh_intp%intp(destination_time, self%slh_raw)
    ENDIF

    IF (self%l_osc_exist) THEN
      ! get new uosc
      CALL self%uosc_intp%intp(destination_time, self%uosc_raw)
      ! get new vosc
      CALL self%vosc_intp%intp(destination_time, self%vosc_raw)
    END IF

  END SUBROUTINE read_wave_forcing__get_new_rawdata


  !>
  !! The forcing state vector is updated
  !!
  !! The forcing state vector is updated by copying the raw data fields to the
  !! corresponding state vars. Additional diagnostic fields are updated on the basis
  !! of the recent raw data fields, which results in a fully updated forcing state vector.
  !!
  SUBROUTINE read_wave_forcing__update_forcing (self, destination_time, u10m, v10m, &
    &                                           sp10m, dir10m, sic, slh, uosc, vosc, &
    &                                           sp_osc, dir_osc, ice_free_mask_c)

    CHARACTER(len=*), PARAMETER :: routine = modname//':read_wave_forcing__update_forcing'

    CLASS(t_read_wave_forcing), TARGET     :: self
    TYPE(datetime), POINTER, INTENT(IN)    :: destination_time ! validity time of the new data
    REAL(wp),                INTENT(INOUT) :: u10m(:,:), v10m(:,:)     ! zonal and meridional wind
                                                                       ! components at 10m asl
    REAL(wp),                INTENT(INOUT) :: sp10m(:,:)               ! wind speed at 10m asl
    REAL(wp),                INTENT(INOUT) :: dir10m(:,:)              ! wind direction at 10m asl
    REAL(wp),                INTENT(INOUT) :: sic(:,:)                 ! sea ice fraction
    REAL(wp),                INTENT(INOUT) :: slh(:,:)                 ! sea level height
    REAL(wp),                INTENT(INOUT) :: uosc(:,:), vosc(:,:)     ! ocean surface currents
    REAL(wp),                INTENT(INOUT) :: sp_osc(:,:)              ! ocean surface current velocity
    REAL(wp),                INTENT(INOUT) :: dir_osc(:,:)             ! ocean surface current direction
    INTEGER,                 INTENT(INOUT) :: ice_free_mask_c(:,:)     ! ice mask

    ! get new forcing data (read from file)
    CALL self%get_new_rawdata(destination_time)

    ! check if the size of the raw data field matches with the
    ! corresponding field in the forcing state
    IF (self%l_wind_exist) THEN
      CALL check_matching_size(arr1=self%u10m_raw, arr2=u10m)
      CALL check_matching_size(arr1=self%v10m_raw, arr2=v10m)
    END IF

    IF (self%l_ice_exist) THEN
      CALL check_matching_size(arr1=self%sic_raw,  arr2=sic)
    ENDIF

    IF (self%l_slh_exist) THEN
      CALL check_matching_size(arr1=self%slh_raw,  arr2=slh)
    ENDIF

    IF (self%l_osc_exist) THEN
      CALL check_matching_size(arr1=self%uosc_raw, arr2=uosc)
      CALL check_matching_size(arr1=self%vosc_raw, arr2=vosc)
    END IF

    ! copy ray data fields to the forcing state vector
    !
!$OMP PARALLEL
    IF (self%l_wind_exist) THEN
      CALL copy(src=self%u10m_raw(:,1,:,1), dest=u10m)
      CALL copy(src=self%v10m_raw(:,1,:,1), dest=v10m)
    END IF

    IF (self%l_ice_exist) THEN
      CALL copy(src=self%sic_raw  (:,1,:,1), dest=sic)
    ENDIF

    IF (self%l_slh_exist) THEN
      CALL copy(src=self%slh_raw  (:,1,:,1), dest=slh)
    ENDIF

    IF (self%l_osc_exist) THEN
      CALL copy(src=self%uosc_raw (:,1,:,1), dest=uosc)
      CALL copy(src=self%vosc_raw (:,1,:,1), dest=vosc)
    END IF
!$OMP END PARALLEL

    ! update additional diagnostic fields
    !

    ! update wind speed and direction
    CALL update_speed_and_direction(p_patch = self%p_patch, &  ! IN
      &                               u     = u10m,         &  ! IN
      &                               v     = v10m,         &  ! IN
      &                              sp     = sp10m,        &  ! OUT
      &                              dir    = dir10m)          ! OUT

    ! update ocean current velocity and direction
    CALL update_speed_and_direction(p_patch = self%p_patch, &  ! IN
      &                               u     = uosc,         &  ! IN
      &                               v     = vosc,         &  ! IN
      &                              sp     = sp_osc,       &  ! OUT
      &                              dir    = dir_osc)         ! OUT

    ! update ice-free mask
    CALL update_ice_free_mask(p_patch       = self%p_patch,   & ! IN
      &                       sea_ice_c     = sic,            & ! IN
      &                       ice_free_mask = ice_free_mask_c)  ! OUT

  CONTAINS

    SUBROUTINE check_matching_size (arr1, arr2)
      REAL(wp), INTENT(IN) :: arr1(:,:,:,:)
      REAL(wp), INTENT(IN) :: arr2(:,:)

      CHARACTER(len=*), PARAMETER :: routine = modname//':check_matching_size'

      IF (SIZE(arr1,1) /= SIZE(arr2,1) .OR. SIZE(arr1,3) /= SIZE(arr2,2)) THEN
        WRITE(message_text,'(a)') 'SIZE mismatch between source and destination array'
        CALL finish(routine, message_text)
      ENDIF
    END SUBROUTINE
  END SUBROUTINE read_wave_forcing__update_forcing


  !>
  !! Destruct reader of wave forcing data
  !!
  !! Destruct the reader object which reads in wave forcing data from file
  !!
  SUBROUTINE read_wave_forcing__deinit (self)

  !  CHARACTER(len=*), PARAMETER :: routine = modname//':read_wave_forcing__deinit'

    CLASS(t_read_wave_forcing)              :: self

    self%p_patch => NULL()

    ! Destruct reader
    IF (self%l_wind_exist) CALL self%uv10_reader%deinit()
    IF (self%l_ice_exist)  CALL self%sic_reader%deinit()
    IF (self%l_slh_exist)  CALL self%slh_reader%deinit()
    IF (self%l_osc_exist)  CALL self%osc_reader%deinit()

    CALL DO_DEALLOCATE(self%u10m_raw)
    CALL DO_DEALLOCATE(self%v10m_raw)
    CALL DO_DEALLOCATE(self%sic_raw)
    CALL DO_DEALLOCATE(self%slh_raw)
    CALL DO_DEALLOCATE(self%uosc_raw)

  END SUBROUTINE read_wave_forcing__deinit

END MODULE mo_wave_forcing
