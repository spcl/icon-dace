!
! mo_art_external_init_radioact
! This module initialises the emission routine of radioactive particles.
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

MODULE mo_art_external_init_radioact
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_io_units,                      ONLY: find_next_free_unit
  USE mo_exception,                     ONLY: finish
  USE mtime,                            ONLY: datetime, timedelta, max_datetime_str_len,        &
                                          &   newTimedelta,newDatetime, getPTStringFromSeconds, &
                                          &   datetimeToString, deallocateDatetime,             &
                                          &   deallocateTimedelta, OPERATOR(+)
  USE, INTRINSIC :: iso_c_binding,      ONLY: c_int64_t
  USE mo_impl_constants,                ONLY: SUCCESS
  USE mo_key_value_store,               ONLY: t_key_value_store
! ART
  USE mo_art_impl_constants,            ONLY: IART_VARNAMELEN, std_radioact_names
  USE mo_art_pntSrc_types,              ONLY: t_art_pntSrc

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_external_init_radioact'

  PUBLIC :: art_extinit_radioact_pntSrc
  PUBLIC :: art_get_nsource_radio_fromfile

  CONTAINS

SUBROUTINE art_get_nsource_radio_fromfile(radioactfile, nsources)
!<
! SUBROUTINE art_get_nsource_radio_fromfile
! This subroutine gets the number of active source scenarios from an
! external input file specified via cart_radioact_file.
! Based on: -
! Part of Module: mo_art_external_init_radioact
! Author: Daniel Rieger, KIT
! Initial Release: 2017-04-07
! Modifications:
! YYYY-MM-DD: <name>, <instituttion>
! - ...
!>
  CHARACTER(LEN=*),INTENT(in) :: &
    &  radioactfile                !< Absolute path + filename of input file for radioactive emissions
  INTEGER, INTENT(out)        :: &
    &  nsources                    !< Number of source scenarios in input file
! Local Variables
  CHARACTER(LEN=40)           :: &
    &  textdummy                   !< Dummy to read unneeded text from input file
  INTEGER                     :: &
    &  iunit,                    & !< Filenumber of the opened dataset for read instructions
    &  iostat,                   & !< File io status
    &  nlines                      !< Loop counter

  nsources = 0

  iunit = find_next_free_unit(100,1000)
  IF (iunit < 0) THEN  
    CALL finish(TRIM(routine)//':art_get_nsource_radio_fromfile',   &
      &         'Failed call to find_next_free_unit.')
  END IF

  OPEN( iunit, file=TRIM(radioactfile), FORM='FORMATTED',  &
    &   STATUS='OLD', ACTION='READ', IOSTAT=iostat )
  IF (iostat /= 0) THEN
    CALL finish(TRIM(routine)//':art_get_nsource_radio_fromfile',   &
      &         'ART: Radioactive input file not found: '//TRIM(radioactfile))
  ELSE
    ! Read the meta information header of the file
    DO nlines = 1, 7
      IF (nlines .EQ. 7) THEN !< Line 7 in the file gives the number active source scenarios
        READ (iunit,'(A40,I10)') textdummy, nsources
      ELSE
        READ(iunit,*)
      ENDIF
    ENDDO

    CLOSE(iunit)
  ENDIF

END SUBROUTINE art_get_nsource_radio_fromfile
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_extinit_radioact_pntSrc(p_patch, z_ifc, tc_dt_model, tc_exp_refdate, radioactfile, dict_tracer, &
  &                                    pntSrc_container, nsources_offset, nsources_radio)
!<
! SUBROUTINE art_extinit_radioact_pntSrc
! This subroutine initialises the emission routine of radioactive particles
! Based on: -
! Part of Module: mo_art_external_init_radioact
! Author: Daniel Rieger, KIT
! Initial Release: 2013-02-28
! Modifications:
! YYYY-MM-DD: <name>, <instituttion>
! - ...
!>
  TYPE(t_patch), INTENT(in)            :: &
    &  p_patch                              !< Current domain
  REAL(wp), INTENT(in)                 :: &
    &  z_ifc(:,:,:)                         !< height of half levels
  TYPE(timedelta), POINTER, INTENT(in) :: &
    &  tc_dt_model                          !< Model timestep
  TYPE(datetime), POINTER, INTENT(in)  :: &
    &  tc_exp_refdate                       !< Experiment reference date
  CHARACTER(LEN=*),INTENT(in)          :: &
    &  radioactfile                         !< Absolute path + filename of input file for radioactive emissions
  TYPE(t_key_value_store), INTENT(in)  :: &
    &  dict_tracer                          !< Tracer index dictionary
  TYPE(t_art_pntSrc), INTENT(inout)    :: &
    &  pntSrc_container(:)                  !< Container with all point sources
  INTEGER, INTENT(in)                  :: &
    &  nsources_offset,                   & !< Sources defined by XML file (here only offset in container)
    &  nsources_radio                       !< Sources to be filled by input file
! Local variables
  TYPE(datetime), POINTER              :: &
    &  startTime, endTime                   !< Start and end time of source scenario
  TYPE(timedelta), POINTER             :: &
    &  StartDeltaFromRefDate,             & !< Timedelta between reference date and start of release
    &  EndDeltaFromRefDate                  !< Timedelta between reference date and end of release
  CHARACTER(LEN=max_datetime_str_len)  :: &
    &  cstartTime, cendTime,               & !< Start and end time of source scenario, format: 'YYYY-MM-DDTHH:MM:SS'
    &  PTStringSecondsAfterRefDate           !< Time after reference date in PT string format
  REAL(wp)                             :: &
    &  src_lon, src_lat,                  & !< Longitude and latitude read from file
    &  src_height,                        & !< Height of source
    &  std_src_strngth,                   & !< standard source strength
    &  src_strngth(9)                       !< source strength for each of the 9 species
  CHARACTER(LEN=10)                    :: &
    &  c_tstart,                          & !< Start time of the emission (from initialization time)
    &  c_tend                               !< End time of the emission (from initialization time)
  CHARACTER(LEN=40)                    :: &
    &  textdummy                            !< Dummy to read unneeded text from input file
  CHARACTER(LEN=IART_VARNAMELEN)       :: &
    &  src_id                               !< Name of source
  INTEGER                              :: &
    &  iunit,                             & !< Filenumber of the opened dataset for read instructions
    &  iostat, ierror,                    & !< File io status, error return value
    &  nlines, jsource, n_isotop,         & !< Loop counter
    &  isource,                           & !< Internal source counter
    &  isec, imin, ihour,                 & !< start/end of emission (seconds/minutes/hours)
    &  idx_tracer(9)                        !< Index of radioactive species in tracer container

  DO n_isotop = 1, 9
    CALL dict_tracer%get(TRIM(ADJUSTL(std_radioact_names(n_isotop))),idx_tracer(n_isotop),ierror)
      IF (ierror /= SUCCESS) CALL finish (TRIM(routine)//'art_extinit_radioact_pntSrc', &
                               &          TRIM(ADJUSTL(std_radioact_names(n_isotop)))//' not found in dictionary.')
  ENDDO

  iunit = find_next_free_unit(100,1000)
  IF (iunit < 0) THEN  
    CALL finish(TRIM(routine)//':art_extinit_radioact_pntSrc',   &
      &         'Failed call to find_next_free_unit.')
  ENDIF

  OPEN( iunit, file=TRIM(radioactfile), FORM='FORMATTED',  &
    &   STATUS='OLD', ACTION='READ', IOSTAT=iostat )
  IF (iostat /= 0) THEN
    CALL finish(TRIM(routine)//':art_extinit_radioact_pntSrc',   &
      &         'ART: Radioactive input file not found: '//TRIM(radioactfile))
  ELSE
    ! Discard the meta information header of the file
    DO nlines = 1, 8
        READ(iunit,*)
    ENDDO

    ! Line 9 contains the standard source strength that is needed:
    READ (iunit,'(A40,F20.0)') textdummy, std_src_strngth

    ! Line 10:
    READ(iunit,*)

    ! Loop over the number of active source scenarios
    DO jsource = 1, nsources_radio

      READ(iunit,'(6X,F10.5,6X,F10.5,5X,F11.5,6X,A10,6X,A10,11X,A15)')               &
        &  src_lon, src_lat, src_height, c_tstart, c_tend, src_id
      READ(iunit,'(11X,15(1X,E9.3))') (src_strngth(n_isotop),n_isotop=1,9)

      READ(c_tstart(9:10),'(I2)') isec
      READ(c_tstart(6:7), '(I2)') imin
      READ(c_tstart(1:4), '(I4)') ihour
      ! Total number of seconds after begin
      isec = ihour*3600 + imin*60 + isec
      CALL getPTStringFromSeconds(INT(isec,c_int64_t),PTStringSecondsAfterRefDate)
      StartDeltaFromRefDate => newTimedelta(TRIM(PTStringSecondsAfterRefDate))
      startTime => newDatetime(tc_exp_refdate, errno=ierror)
      IF(ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_extinit_radioact_pntSrc', &
                              &         'Could not create datetime object from tc_exp_refdate')
      ! Add delta from file to start date of simulation
      startTime = startTime + StartDeltaFromRefDate
      CALL datetimeToString(startTime, cstartTime)
      CALL deallocateDatetime(startTime)
      CALL deallocateTimedelta(StartDeltaFromRefDate)

      READ(c_tend(9:10),  '(I2)') isec
      READ(c_tend(6:7),   '(I2)') imin
      READ(c_tend(1:4),   '(I4)') ihour
      ! Total number of seconds after begin
      isec = ihour*3600 + imin*60 + isec
      CALL getPTStringFromSeconds(INT(isec,c_int64_t),PTStringSecondsAfterRefDate)
      EndDeltaFromRefDate => newTimedelta(TRIM(PTStringSecondsAfterRefDate))
      endTime => newDatetime(tc_exp_refdate, errno=ierror)
      IF(ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_extinit_radioact_pntSrc', &
                              &         'Could not create datetime object from tc_exp_refdate')
      ! Add delta from file to start date of simulation
      endTime = endTime + EndDeltaFromRefDate
      CALL datetimeToString(endTime, cendTime)
      CALL deallocateDatetime(endTime)
      CALL deallocateTimedelta(EndDeltaFromRefDate)

      ! Add a source for each isotope
      DO n_isotop = 1, 9
        isource = nsources_offset + (jsource-1)*9 + n_isotop
        src_strngth(n_isotop) = src_strngth(n_isotop) * std_src_strngth / 3600._wp ! Bq h-1 -> Bq s-1

        CALL pntSrc_container(isource)%init(tc_dt_model, tc_exp_refdate, p_patch, z_ifc,                             &
          &                                 TRIM(ADJUSTL(src_id))//'-'//TRIM(ADJUSTL(std_radioact_names(n_isotop))), &
          &                                 src_lon, src_lat, src_height,                                            &
          &                                 idx_tracer(n_isotop), src_strngth(n_isotop),                             &
          &                                 cstartTime, cendTime)

        CALL pntSrc_container(isource)%print
      ENDDO !n_isotop
    ENDDO !jsource

    CLOSE(iunit)
  ENDIF

END SUBROUTINE art_extinit_radioact_pntSrc
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_external_init_radioact
