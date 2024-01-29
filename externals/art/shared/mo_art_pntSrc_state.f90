!
! mo_art_pntSrc_state
! This module initializes the data structure required by the pntSrc component
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

MODULE mo_art_pntSrc_state
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_key_value_store,               ONLY: t_key_value_store
  USE mo_exception,                     ONLY: message, message_text, finish
  USE mo_impl_constants,                ONLY: SUCCESS
  USE mo_math_constants,                ONLY: pi
  USE mtime,                            ONLY: timedelta, datetime
  USE mo_tracer_metadata_types,         ONLY: t_aero_meta, t_chem_meta
  USE mo_var_list,                      ONLY: t_var_list_ptr
  USE mo_var_metadata_types,            ONLY: t_var_metadata_dynamic
! ART
  USE mo_art_data,                      ONLY: t_art_data
  USE mo_art_modes_linked_list,         ONLY: p_mode_state, t_mode
  USE mo_art_modes,                     ONLY: t_fields_2mom
  USE mo_art_read_xml,                  ONLY: t_xml_file, art_open_xml_file, art_close_xml_file,  &
                                          &   art_get_childnumber_xml, art_read_elements_xml
  USE mo_art_external_init_radioact,    ONLY: art_get_nsource_radio_fromfile, art_extinit_radioact_pntSrc
  USE mo_art_create_filenames,          ONLY: art_check_filename
  USE mo_art_impl_constants,            ONLY: UNDEF_REAL_ART, UNDEF_INT_ART
  USE mo_art_string_tools,              ONLY: key_value_storage_as_string

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_pntSrc_state'

  PUBLIC :: art_emiss_init_pntSrc

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_emiss_init_pntSrc(p_patch, z_ifc, tc_dt_model, tc_exp_refdate,                     &
  &                              cpntSrc_xml_file, iart_radioact, radioactfile, p_prog_list,      &
  &                              p_art_data)
!<
! SUBROUTINE art_emiss_init_pntSrc
! This subroutine initializes the data structure required by mo_art_emission_pntSrc
! Based on: Werchner (2016) - Bachelorthesis, KIT
! Part of Module: mo_art_pntSrc_state
! Author: Sven Werchner, Daniel Rieger, KIT
! Initial Release: 2017-01-25
! Modifications:
! 2018-08-28: Lukas Muser, KIT
! - read emission profile from xml
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  TYPE(t_patch), INTENT(in)      :: &
    &  p_patch                        !< Current domain
  REAL(wp), INTENT(in)           :: &
    &  z_ifc(:,:,:)                   !< height of half levels
  TYPE(timedelta), POINTER, INTENT(in) :: &
    &  tc_dt_model                    !< Model timestep
  TYPE(datetime), POINTER, INTENT(in)  :: &
    &  tc_exp_refdate                 !< Experiment reference date
  CHARACTER(LEN=*),INTENT(in)    :: &
    &  cpntSrc_xml_file               !< XML file containing point source information
  INTEGER, INTENT(in)            :: &
    &  iart_radioact                  !< Radioactive emission configuration
  TYPE(t_var_list_ptr),INTENT(in)    :: &
    &  p_prog_list                    !< current list: prognostic
  CHARACTER(LEN=*),INTENT(in)    :: &
    &  radioactfile                   !< File with radioactive emission information (official format)
  TYPE(t_art_data),INTENT(inout) :: &
    &  p_art_data                     !< ART data container
! Local variables
  TYPE(t_key_value_store)        :: &
    &  key_value_store                !< Storage container for point source metadata from XML
  TYPE(t_xml_file)               :: &
    &  pntSrc_xmlfile                 !< point source XML file
  REAL(wp)                       :: &
    &  lon,                         & !< longitude of the emission source
    &  lat,                         & !< latitude of the emission source
    &  height,                      & !< height of the emission source above ground (m)
    &  height_bot,                  & !< bottom height of the emission source above ground (m)
    &  source_strength,             & !< Source strength of tracer per s-1
    &  emiss_rate0,                 & !< for aerosol tracer (modal distribution) emission -> emission rate of number concentration
    &  dg3_emiss, sigma_emiss         !< median diameter and std of emitted tracer mass ratio
  INTEGER                        :: &
    &  itr,                         & !< Index of tracer to apply this source scenario
    &  itr0,                        & !< index of corresponding nmb_conc for aero tracer
    &  jsource,                     & !< Counter for sources
    &  nsources,                    & !< Total number of point sources
    &  nsources_radio,              & !< Total number of point sources from radioactive input file
    &  ierror
  CHARACTER(:), ALLOCATABLE :: &
    &  id,                          & !< Name of source
    &  tracer_name,                 & !< Name of tracer the source corresponds to
    &  unit                           !< unit of emission flux
  CHARACTER(LEN=3)               :: &
    &  jsource_str                    !< jsource as character
  CHARACTER(:), ALLOCATABLE  :: &
    &  startTime, endTime,          & !< Start and end time of source scenario, format: 'YYYY-MM-DDTHH:MM:SS'
    &  emiss_profile                  !< arithmetic expression of emission profile
  LOGICAL                        :: &
    &  lexist

  CALL art_check_filename(TRIM(cpntSrc_xml_file), lrequired=.TRUE.)

  CALL art_open_xml_file(TRIM(cpntSrc_xml_file),pntSrc_xmlfile)

  ! read xml file with sources
  CALL art_get_childnumber_xml(pntSrc_xmlfile, '/sources/', nsources)
  WRITE (message_text,'(A31,I4)') 'ART: Number of active sources: ',nsources
  CALL message (TRIM(routine)//':art_emiss_init_pntSrc', message_text)
  p_art_data%pntSrc%nsources = nsources
  
  ! check radioactive input file for sources
  IF ( iart_radioact == 1 .AND. TRIM(radioactfile) /= '') THEN
    CALL art_check_filename(TRIM(radioactfile), lrequired=.FALSE., lexist=lexist)
    IF (lexist) THEN
      CALL art_get_nsource_radio_fromfile(TRIM(radioactfile), nsources_radio)
      WRITE (message_text,'(A47,I4)') 'ART: Number of additional radioactive sources: 9 x ',nsources_radio
      CALL message (TRIM(routine)//':art_emiss_init_pntSrc', message_text)
      p_art_data%pntSrc%nsources = nsources + (nsources_radio * 9) ! each source is active for 9 tracer
    ELSE
      WRITE (message_text,'(A47,I4)') 'ART: WARNING: cart_radioact_file = '//TRIM(radioactfile)//' not found.'
      CALL message (TRIM(routine)//':art_emiss_init_pntSrc', message_text)
    ENDIF
  ENDIF

  ! create point source container
  ALLOCATE(p_art_data%pntSrc%p(p_art_data%pntSrc%nsources))

  DO jsource = 1, nsources

    CALL key_value_store%init(.FALSE.)

    ! Read name and metadata of current tracer in XML file
    WRITE(jsource_str,'(I3)') jsource
    CALL art_read_elements_xml(pntSrc_xmlfile,'/sources/*['//TRIM(ADJUSTL(jsource_str))//']/',key_value_store)

    ! get information of current source -> First step: No default values
    CALL key_value_storage_as_string(key_value_store,'name',id,ierror)
      IF(ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_emiss_init_pntSrc', &
                              &         'id not found in XML file.')
    CALL key_value_store%get('lon',lon,ierror)
      IF(ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_emiss_init_pntSrc', &
                              &         'lon not found in XML file.')
    CALL key_value_store%get('lat',lat,ierror)
      IF(ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_emiss_init_pntSrc', &
                              &         'lat not found in XML file.')
    CALL key_value_storage_as_string(key_value_store,'substance',tracer_name,ierror)
      IF(ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_emiss_init_pntSrc', &
                              &         'substance not found in XML file.')
    CALL key_value_store%get('source_strength',source_strength,ierror)
      IF(ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_emiss_init_pntSrc', &
                              &         'source_strength not found in XML file.')
    CALL key_value_store%get('height',height,ierror)
      IF(ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_emiss_init_pntSrc', &
                              &         'height not found in XML file.')

    CALL key_value_store%get('height_bot',height_bot,ierror)
      IF (ierror /= SUCCESS) height_bot = UNDEF_REAL_ART

    CALL key_value_storage_as_string(key_value_store,'emiss_profile',emiss_profile,ierror)
      IF (ierror /= SUCCESS) emiss_profile = ''

    CALL key_value_store%get('dg3_emiss',dg3_emiss,ierror)
      IF(ierror /= SUCCESS) dg3_emiss = UNDEF_REAL_ART
    CALL key_value_store%get('sigma_emiss',sigma_emiss,ierror)
      IF(ierror /= SUCCESS) sigma_emiss = UNDEF_REAL_ART

    CALL key_value_storage_as_string(key_value_store,'unit',unit,ierror)
      IF(ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_emiss_init_pntSrc', &
                              &         'unit not found in XML file.')
    CALL convert_unit_pntSrc(TRIM(unit),source_strength)

    CALL p_art_data%dict_tracer%get(TRIM(tracer_name),itr,ierror)
      IF(ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_emiss_init_pntSrc', &
                              &         'Tracer '//TRIM(tracer_name)//' not found in dictionary.')

    CALL tracer_metadata_for_pntSrc(p_prog_list, TRIM(tracer_name), itr, dg3_emiss, sigma_emiss,  &
      &                             source_strength, itr0, emiss_rate0)


    CALL key_value_storage_as_string(key_value_store,'startTime',startTime,ierror)
      IF (ierror /= SUCCESS) THEN
        startTime = '1582-10-15T00:00:00'
        WRITE (message_text,*) 'ART: WARNING: startTime of  '//TRIM(id)//' not found: Setting to 1582-10-15T00:00:00.'
        CALL message (TRIM(routine)//':art_emiss_init_pntSrc', message_text)
      ENDIF
    CALL key_value_storage_as_string(key_value_store,'endTime',endTime,ierror)
      IF (ierror /= SUCCESS) THEN
        endTime = '9999-12-31T00:00:00'
        WRITE (message_text,*) 'ART: WARNING: endTime of  '//TRIM(id)//' not found: Setting to 9999-12-31T00:00:00.'
        CALL message (TRIM(routine)//':art_emiss_init_pntSrc', message_text)
      ENDIF

    CALL p_art_data%pntSrc%p(jsource)%init(tc_dt_model, tc_exp_refdate, p_patch, z_ifc, id,       &
      &                                    lon, lat, height, itr, source_strength,                &
      &                                    startTime, endTime,                                    &
      &                                    emiss_profile=TRIM(ADJUSTL(emiss_profile)),            &
      &                                    height_bot=height_bot, itr0=itr0,                      &
      &                                    emiss_rate0=emiss_rate0 )


    CALL p_art_data%pntSrc%p(jsource)%print
    WRITE (message_text,'(A16,A)')     'emiss_profile:  ',TRIM(ADJUSTL(emiss_profile))
    CALL message ('', message_text)

    CALL key_value_store%destruct
  ENDDO

  CALL art_close_xml_file(pntSrc_xmlfile)

  IF ( iart_radioact == 1 .AND. TRIM(radioactfile) /= '') THEN
    CALL art_extinit_radioact_pntSrc(p_patch, z_ifc, tc_dt_model, tc_exp_refdate, TRIM(radioactfile),  &
      &                              p_art_data%dict_tracer, p_art_data%pntSrc%p, nsources, nsources_radio)
  ENDIF

END SUBROUTINE art_emiss_init_pntSrc
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE convert_unit_pntSrc(unit,source_strength)
!<
! SUBROUTINE convert_unit_pntSrc
! This subroutine converts the source strength from a given unit to kg s-1
! Based on: -
! Part of Module: mo_art_pntSrc_state
! Author: Daniel Rieger, KIT
! Initial Release: 2017-04-11
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  CHARACTER(LEN=*), INTENT(in)   :: &
    &  unit                           !< unit of emission flux
  REAL(wp), INTENT(inout)        :: &
    &  source_strength                !< Source strength of tracer in: unit , out: s-1

! Maybe we could check here also if the unit corresponds to the tracer unit. But this takes some 
! effort as we first need to get the tracer entry from the container.

  SELECT CASE(TRIM(unit))
    CASE('kg s-1','Bq s-1','s-1')
      ! Nothing to do
    CASE('kg h-1','Bq h-1','h-1')
      source_strength = source_strength/3600._wp
    CASE('kg d-1','Bq d-1','d-1')
      source_strength = source_strength/86400._wp
    CASE DEFAULT
      CALL finish(TRIM(routine)//':convert_unit_pntSrc', &
        &         'Unit '//TRIM(unit)//' is not supported for point sources.')
  END SELECT

END SUBROUTINE convert_unit_pntSrc
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE tracer_metadata_for_pntSrc(p_prog_list,tr_name,itr,dg3,sigma,emiss_rate,  &
  &                                   itr0,emiss_rate0)
!<
! SUBROUTINE tracer_metadata_for_pntSrc
! This subroutine gets information about the emitted tracer's metadata
! Based on: -
! Part of Module: mo_art_pntSrc_state
! Author: Lukas Muser, KIT
! Initial Release: 2018-09-12
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  TYPE(t_var_list_ptr), INTENT(in)    :: &
    &  p_prog_list                     !< list of prognostic variables
  CHARACTER(LEN=*), INTENT(in)    :: &
    &  tr_name                         !< name of tracer that is emitted by pntSrc
  INTEGER,          INTENT(in)    :: &
    &  itr                             !< index of (mass-) tracer in tracer container
  REAL(wp),         INTENT(in)    :: &
    &  dg3, sigma                      !< for aerosol tracer -> median diameter and standard
                                       !  deviation of emitted tracer distribution
  REAL(wp),         INTENT(inout) :: &
    &  emiss_rate                      !< emission rate of pntSrc in tracer specific units
  INTEGER,          INTENT(out)   :: & 
    &  itr0                            !< index of corresponding nmb_conc for aero tracer
  REAL(wp),         INTENT(out)   :: &
    &  emiss_rate0                     !< for aerosol tracer (modal distribution) emission 
                                       !  -> emission rate of number concentration
! Local variables
  TYPE(t_var_metadata_dynamic),POINTER     :: &
    & info_dyn                         !< returns reference to dynamic tracer metadata
  TYPE(t_mode),     POINTER       :: &
    &  current_mode                    !< pointer to loop through mode structure
  REAL(wp)                        :: &
    &  unit_conv,                    & !< factor that converts aerosol tracer emission into amug/kg
    &  rho                             !< density of aerosol tracer
  INTEGER                         :: &
    &  ierror, iv
  LOGICAL                         :: &
    &  lcalc_2mom

  lcalc_2mom  = .TRUE.
  itr0        = UNDEF_INT_ART   ! Default for all tracers other then aero tracer in mode t_fields_2mom
  emiss_rate0 = UNDEF_REAL_ART  ! Default for all tracers other then aero tracer in mode t_fields_2mom

  !CALL get_tracer_info_dyn_by_idx(p_prog_list, itr, info)
  DO iv = 1, p_prog_list%p%nvars
    IF(p_prog_list%p%vl(iv)%p%info%ncontained /= itr) CYCLE
    info_dyn => p_prog_list%p%vl(iv)%p%info_dyn
  END DO
  
  SELECT TYPE(tracer_info => info_dyn%tracer)
    CLASS IS(t_aero_meta)
      ! Loop through modes and find itr0 -> index of number concentration
      current_mode => p_mode_state(1)%p_mode_list%p%first_mode
      DO WHILE(ASSOCIATED(current_mode))
        ! Select type of mode
        SELECT TYPE (fields=>current_mode%fields)
          CLASS IS (t_fields_2mom)
            IF (TRIM(tracer_info%mode) == TRIM(fields%name)) THEN
              itr0 = fields%itr0
            ENDIF
            IF (lcalc_2mom) THEN
              IF (dg3 == UNDEF_REAL_ART)   CALL finish(TRIM(routine)//':init_pntSrc',                     &
                                        &  'dg3_emiss has to be defined in pntSrc.xml for aero tracer')
              IF (sigma == UNDEF_REAL_ART) CALL finish(TRIM(routine)//':init_pntSrc',                     &
                                        &  'sigma_emiss has to be defined in pntSrc.xml for aero tracer')
              CALL tracer_info%opt_meta%get('rho', rho, ierror)
              IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':tracer_metadata_for_pntSrc',            &
                                        &  'rho not available for tracer '//TRIM(info_dyn%tracer%name)//'.')
              emiss_rate0 = 6.0_wp / pi / rho                             &
                &         * EXP( 4.5_wp * (LOG(sigma)**2.0_wp)) / (dg3**3.0_wp)     &
                &         * emiss_rate
              unit_conv = 1.0e9_wp                           ! convert kg/kg in amug/kg
              lcalc_2mom = .FALSE.
            END IF
          CLASS DEFAULT
            ! no unit conversion for pollen and radiocative species
            ! itr0 is UNDEF_INT_ART
            unit_conv = 1.0_wp
        END SELECT !fields
        current_mode => current_mode%next_mode
      ENDDO !ASSOCIATED(current_mode)
      
      emiss_rate  = unit_conv * emiss_rate

    CLASS IS(t_chem_meta)
      unit_conv   = 1.0_wp
      emiss_rate  = unit_conv * emiss_rate

    CLASS DEFAULT
      unit_conv = 1.0_wp
      emiss_rate  = unit_conv * emiss_rate
      WRITE (message_text,*) 'WARNING: ART: Unkown metadata type of tracer ',tr_name
      CALL message (TRIM(routine)//':tracer_metadata_for_pntSrc', message_text)
  END SELECT

END SUBROUTINE tracer_metadata_for_pntSrc
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_pntSrc_state
