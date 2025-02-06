!
! mo_art_emiss_state
! This module provides subroutines for initialising the
! emission metadata from the xml file given by namelist parameter
!
!
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

MODULE mo_art_emiss_state
  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: finish, message, message_text
  USE mo_impl_constants,        ONLY: SUCCESS
  USE mo_art_config,            ONLY: t_art_config, art_config, IART_PATH_LEN
  USE mo_run_config,            ONLY: ntracer
  USE mo_var_list,              ONLY: t_var_list_ptr
  USE mo_var_metadata_types,    ONLY: t_var_metadata, t_var_metadata_dynamic
  USE mo_tracer_metadata_types, ONLY: t_chem_meta, t_aero_meta
  USE mo_impl_constants,        ONLY: SUCCESS
  USE mo_key_value_store,       ONLY: t_key_value_store
  !EXTERNALS
  USE mtime,                    ONLY: max_datetime_str_len, &
                                  &   datetimeToPosixString, &
                                  &   getNoOfDaysInYearDateTime
  !ART
  USE mo_art_data,              ONLY: p_art_data, t_art_data
  USE mo_art_atmo_data,         ONLY: t_art_atmo
  USE mo_art_prescribed_types,  ONLY: t_art_emiss_prescribed
  USE mo_art_emiss_types,       ONLY: t_art_emiss_type_container, &
                                  &   t_art_emiss2tracer
  USE mo_art_io_constants,      ONLY: art_get_all_emission_elements
  USE mo_art_create_filenames,  ONLY: art_create_filenames
  USE mo_art_read_xml,          ONLY: t_xml_file,                  &
                                  &   art_read_emission_dataset,   &
                                  &   art_read_elements_xml,       &
                                  &   art_get_childnumber_xml,     &
                                  &   art_open_xml_file,           &
                                  &   art_close_xml_file
  USE mo_art_bvoc,              ONLY: pftn, sel_bccnum
  USE mo_art_read_extdata,      ONLY: art_read_PFTdata
  USE mo_art_io_constants,      ONLY: prescr_emiss_elements,   &
                                  &   iemiss, IEMISS_BIO,      &
                                  &   STORAGE_ATTR_SEP
  USE mo_art_impl_constants,    ONLY: IART_VARNAMELEN, IART_EMISS2TRACER

  USE mo_art_modes_linked_list, ONLY: t_mode_state, t_mode, new_mode_list, append_mode, t_mode_list
  USE mo_art_modes,             ONLY: t_fields_2mom
  USE mo_art_string_tools,      ONLY: key_value_storage_as_string

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: routine = 'mo_art_emiss_state'

  PUBLIC :: art_init_emissions_from_xml
  PUBLIC :: art_aerosol_assign_emission2tracer

CONTAINS

!
!------------------------------------------------------------------------------------
!

SUBROUTINE art_init_emissions_from_xml(jg,p_prog_list,cart_emiss_xml)
!<
! SUBROUTINE art_init_emissions_from_xml
! This subroutine works as wrapper between mo_art_init and the emission
! initialisation so that the xml files do not have to be opened and closed in
! mo_art_init
! Part of Module: mo_art_emiss_state
! Author: Michael Weimer, KIT
! Initial Release: 2017-07-18
! Modifications:
! yyyy-mm-dd:
! - 
!>
  IMPLICIT NONE
  INTEGER, INTENT(IN)       ::  &
    &  jg              !< patch on which computation is performed
  TYPE(t_var_list_ptr), INTENT(IN) :: &
    &  p_prog_list     !< list of prognostic variables
  CHARACTER(LEN=*), INTENT(IN) :: &
    &  cart_emiss_xml  !< path and file name of the XML file containing the meta
                       !  information of the emissions
  ! local variables
  TYPE(t_xml_file) :: &
    &  tixi_file_em    !< path and handle of the xml file

  IF (cart_emiss_xml /= '') THEN
    CALL art_open_xml_file(TRIM(cart_emiss_xml),tixi_file_em)
  ELSE
    tixi_file_em%path = ""
    tixi_file_em%handle = -1
  END IF

  CALL art_init_tracer_emissions(jg, p_prog_list, &
                    &            tixi_file_em)


  IF (cart_emiss_xml /= '') THEN
    CALL art_close_xml_file(tixi_file_em)
  END IF

END SUBROUTINE art_init_emissions_from_xml
!
!------------------------------------------------------------------------------------
!

SUBROUTINE art_init_tracer_emissions(jg,p_prog_list,tixi_file_emiss)
!<
! SUBROUTINE art_init_tracer_emissions
! This subroutine takes the storage emission elements and stores them in the
! internal emission structure
! Part of Module: mo_art_emiss_state
! Author: Michael Weimer, KIT
! Initial Release: 2016-10-26
! Modifications:
! 2017-07-20: Michael Weimer, KIT
! - included storage as meta data container
! yyyy-mm-dd:
! - 
!>
  IMPLICIT NONE
  INTEGER, INTENT(IN)        ::  &
    &  jg                             !< patch on which computation is performed
  TYPE(t_var_list_ptr), INTENT(IN) :: &
    &  p_prog_list                    !< list of prognostic variables
  TYPE(t_xml_file), INTENT(IN) :: &
    &  tixi_file_emiss                !< emission datasets XML file
  !local variables
  CHARACTER(:), ALLOCATABLE :: &
    &  tracer_name,     &             !< tracer name
    &  dummy_char                     !< dummy string when checking for the
                                      !  existence of emissions in the storage
  INTEGER ::            &
    &  tracer_idx                     !< index of the tracer
  INTEGER ::            &
    &  ierror,          &             !< error when reading meta data
    &  itype,           &             !< type index of emission
    &  iv,              &             !< loop index (list elements)
    &  icur_dataset                   !< index of the prescribed dataset
  TYPE(t_art_emiss_type_container) ::  &
    &  emiss                          !< local element of emission metadata to be added 
                                      !  to the emission storage in p_art_data(jg)
  TYPE(t_key_value_store), POINTER ::  &
    &  meta_storage                   !< pointer to the tracer storage
  CHARACTER(LEN = 3) :: &
    &  sel_bccnum_str                 !< string of sel_bccnum
  CHARACTER(LEN = 5) :: &
    &  nlev_str                       !< string of p_patch%nlev
  CHARACTER(LEN=*), PARAMETER :: &
    & routine = 'mo_art_emiss_state:art_init_tracer_emissions'
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo

  TYPE(t_chem_meta), POINTER :: chem_meta_ptr
  TYPE(t_aero_meta), POINTER :: aero_meta_ptr

  art_atmo => p_art_data(jg)%atmo
  NULLIFY(p_art_data(jg)%ext%land%pft)


  DO iv = 1, p_prog_list%p%nvars

    IF (p_prog_list%p%vl(iv)%p%info_dyn%tracer%lis_tracer) THEN
      tracer_idx = p_prog_list%p%vl(iv)%p%info%ncontained
  
      ! set pointer to the tracer storage or skip the tracers with wrong types
      SELECT TYPE(meta => p_prog_list%p%vl(iv)%p%info_dyn%tracer)
        CLASS IS (t_chem_meta)
          chem_meta_ptr => meta
          meta_storage => chem_meta_ptr%opt_meta
        TYPE IS (t_aero_meta)
          aero_meta_ptr => meta
          meta_storage => aero_meta_ptr%opt_meta
        CLASS DEFAULT
          CYCLE
      END SELECT

      ! get the name of the tracer
      CALL key_value_storage_as_string(meta_storage,'name',tracer_name)
  
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! online biogenic emissions
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CALL meta_storage%get('emiss_onlBIO',emiss%bioonl%idx,ierror)
      IF (ierror /= SUCCESS) THEN
        emiss%bioonl%idx = 0
      ELSE
        ! consistency check of emiss%bioonl%idx
        IF ((emiss%bioonl%idx < 1) .OR. (emiss%bioonl%idx > sel_bccnum)) THEN
          WRITE(sel_bccnum_str,'(A3)') sel_bccnum

          CALL finish(routine, 'Index of online biogenic emissions has to be ' &
                         &   //'between 1 and '//TRIM(sel_bccnum_str)//' for ' &
                         &   //TRIM(tracer_name)//'.' )
        END IF

        ! read inum_levs attribiute if necessary
        IF (.NOT. art_config(jg)%lart_emiss_turbdiff) THEN
          CALL meta_storage%get('inum_levs'//STORAGE_ATTR_SEP   &
                    &         //'emiss_onlBIO',emiss%bioonl%num_emiss_lev,ierror)
          IF (ierror /= SUCCESS)  THEN
            CALL finish(routine,'Attribute inum_levs for element  emiss_onlBIO not found for '  &
                    &         //TRIM(tracer_name))
          END IF

          IF ((emiss%bioonl%num_emiss_lev < 1)   &
             & .OR. (emiss%bioonl%num_emiss_lev > art_atmo%nlev)) THEN
            WRITE(nlev_str,'(I5)') art_atmo%nlev

            CALL finish(routine,'Attribute inum_levs for emiss_onlBIO out of bounds for '  &
                  &           //TRIM(tracer_name)//' (must be between 1 and '  &
                  &           //TRIM(nlev_str)//').')
          END IF
        ELSE
          emiss%bioonl%num_emiss_lev = -1
        END IF

        ! read optional attribute rscaling_factor
        CALL meta_storage%get('rscaling_factor'//STORAGE_ATTR_SEP//'emiss_onlBIO',   &
                 &            emiss%bioonl%scaling_factor,ierror)
        IF (ierror /= SUCCESS) THEN
          emiss%bioonl%scaling_factor = 1._wp
        END IF

        ! read the PFT data (path given in link called "PFT" in
        ! cart_input_folder of art_nml)
        IF (.NOT. ASSOCIATED(p_art_data(jg)%ext%land%pft)) THEN
          ALLOCATE(p_art_data(jg)%ext%land%pft(art_atmo%nproma,art_atmo%nblks,pftn))
          CALL art_read_PFTdata(p_art_data(jg)%ext%land%pft(:,:,:),jg)
        END IF
      
        ALLOCATE(emiss%bioonl%cbio(art_atmo%nproma,art_atmo%nblks))
      
        emiss%bioonl%pft => p_art_data(jg)%ext%land%pft
      END IF
  
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! prescribed emissions
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
      ! determine number of prescribed emissions for this tracer
      emiss%num_types_prescribed = 0

      DO itype = 1,SIZE(prescr_emiss_elements)
        CALL key_value_storage_as_string(meta_storage,prescr_emiss_elements(itype), &
          &                              dummy_char,ierror)
        IF (ierror == SUCCESS) THEN
           IF ((iemiss(itype) == IEMISS_BIO) .AND. (emiss%bioonl%idx > 0)) THEN
               CALL finish(routine,'onlBIO and biogenic emissions at tracer '&
            &                    //TRIM(tracer_name)                         &
            &                    //' are not allowed at the same time.')
          END IF

          emiss%num_types_prescribed = emiss%num_types_prescribed + 1
        END IF
      END DO

  
      IF (emiss%num_types_prescribed > 0) THEN
        IF (TRIM(tixi_file_emiss%path) /= "") THEN
          ALLOCATE(emiss%types(emiss%num_types_prescribed))
          icur_dataset = 1
    
         DO itype = 1,SIZE(prescr_emiss_elements)
            IF(icur_dataset>emiss%num_types_prescribed) CYCLE
            CALL key_value_storage_as_string(meta_storage,prescr_emiss_elements(itype), &
              &                              emiss%types(icur_dataset)%id,ierror)
            IF (ierror == SUCCESS) THEN
              ! read inum_levs attribute if necessary
              IF (.NOT. art_config(jg)%lart_emiss_turbdiff) THEN
                CALL meta_storage%get('inum_levs'//STORAGE_ATTR_SEP  &
                              &     //prescr_emiss_elements(itype),  &
                              &       emiss%types(icur_dataset)%num_emiss_lev,ierror)
                IF (ierror /= SUCCESS) THEN
                  CALL finish(routine,'inum_levs'//STORAGE_ATTR_SEP   &
                         &          //prescr_emiss_elements(itype)  &
                         &          //' not found for '//TRIM(tracer_name))
                END IF

                IF ((emiss%types(icur_dataset)%num_emiss_lev < 1)  &
                  & .OR. (emiss%types(icur_dataset)%num_emiss_lev > art_atmo%nlev)) THEN

                  WRITE(nlev_str,'(I5)') art_atmo%nlev

                  CALL finish(routine,'inum_levs'//STORAGE_ATTR_SEP   &
                        &           //prescr_emiss_elements(itype)  &
                        &           //' out of bounds for '                       &
                        &           //TRIM(tracer_name)//' (must be between 1 and ' &
                        &           //TRIM(nlev_str)//').')
                END IF
              ELSE
                emiss%types(icur_dataset)%num_emiss_lev = -1
              END IF
                
    
              ! read optional rscaling_factor attribute
              CALL meta_storage%get('rscaling_factor'//STORAGE_ATTR_SEP   &
                         &        //prescr_emiss_elements(itype), &
                         &          emiss%types(icur_dataset)%scaling_factor,ierror)
              IF (ierror /= SUCCESS) THEN
                emiss%types(icur_dataset)%scaling_factor = 1._wp
              END IF
    
              emiss%types(icur_dataset)%iType_data = iemiss(itype)
    
              ! read all other metadata from the emission XML file and initialise
              ! the emissions
              CALL art_read_emission_dataset(jg,emiss%types(icur_dataset), &
                                       &     tixi_file_emiss)
    
              CALL art_init_emission_dataset(jg,                        &
                                    &        emiss%types(icur_dataset), &
                                    &        tixi_file_emiss)
    
              icur_dataset = icur_dataset + 1
            END IF

            IF (icur_dataset > emiss%num_types_prescribed) EXIT
          END DO
        ELSE ! path /= ""
          CALL finish(routine,'found prescribed emission in xml files but ' &
                &           //'art_nml:cart_emiss_xml_file has not been given.')
        END IF
      END IF
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! standard emissions
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
      CALL meta_storage%get('emiss_std',emiss%std%val,ierror)
      IF (ierror /= SUCCESS) THEN
        emiss%std%val = 0.0_wp
      ELSE
        IF ((emiss%num_types_prescribed > 0) .OR. (emiss%bioonl%idx > 0)) THEN
          CALL message(routine,'found standard emission but did not consider it ' &
                        &    //'since other emissions were found.')
        ELSE
          ! consistency check of standard value
          IF (emiss%std%val < 0.0_wp) THEN
            CALL finish(routine,'standard emission value has to be greater or '   &
                   &          //'equal zero for '//TRIM(tracer_name))
          END IF
 
          ! read imode attribute
          CALL meta_storage%get('imode'//STORAGE_ATTR_SEP//'emiss_std',   &
                      &         emiss%std%mode,ierror)
          IF (ierror /= SUCCESS) THEN
            CALL finish(routine,'Attribute imode for emiss_std not found for ' &
                  &           //TRIM(tracer_name))
          END IF

          IF (ALL(emiss%std%mode /= (/ 1,2 /)))  THEN
            CALL finish(routine,'Attribute imode for emiss_std has to be 1 or 2 for ' &
                     &        //TRIM(tracer_name))
          END IF
    
          ! read inum_levs attribute if necessary
          IF (.NOT. art_config(jg)%lart_emiss_turbdiff) THEN
            CALL meta_storage%get('inum_levs'//STORAGE_ATTR_SEP//'emiss_std',  &
                      &           emiss%std%num_emiss_lev,ierror)
            IF (ierror /= SUCCESS) THEN
              CALL finish(routine,'Attribute inum_levs for emiss_std'  &
                     &          //'not found for'//TRIM(tracer_name))
            END IF

            IF    ((emiss%std%num_emiss_lev < 1)   &
            & .OR. (emiss%std%num_emiss_lev > art_atmo%nlev)) THEN

              WRITE(nlev_str,'(I5)') art_atmo%nlev

              CALL finish(routine,'Attribute inum_levs for emiss_std out of bounds for ' &
                    &           //TRIM(tracer_name)//' (must be between 1 and ' &
                    &           //TRIM(nlev_str)//').')
            END IF
          ELSE
            emiss%std%num_emiss_lev = -1
          END IF
        END IF
      END IF

      IF ((emiss%num_types_prescribed > 0) .OR. (emiss%bioonl%idx > 0)  &
          &  .OR. (emiss%std%val > 0.0_wp) ) THEN

        ! initialise the emission storage if not already done before

        IF (.NOT. p_art_data(jg)%emiss%is_init) THEN
          CALL p_art_data(jg)%emiss%init(ntracer)
        END IF

        CALL message(routine,'Initialised emissions for tracer ' &
               &           //TRIM(tracer_name))

        CALL p_art_data(jg)%emiss%put(tracer_idx,emiss)

        IF (ALLOCATED(emiss%types)) THEN
          DEALLOCATE(emiss%types)
        END IF

        IF (ALLOCATED(emiss%bioonl%cbio)) THEN
          DEALLOCATE(emiss%bioonl%cbio)
        END IF
      END IF
    END IF ! lis_tracer
       
  END DO

END SUBROUTINE art_init_tracer_emissions
!
!--------------------------------------------------------------------------------------
!
SUBROUTINE art_init_emission_dataset(jg,emiss_t,tixi_file_emiss)
!<
! SUBROUTINE art_init_emission_datasets
! This subroutine checks if the files corresponding to the boundaries of the 
! dataset exist and reads the emissions closest to the current simulation time
! into the structure
! Part of Module: mo_art_emiss_state
! Author: Michael Weimer, KIT
! Initial Release: 2016-10-26
! Modifications:
! yyyy-mm-dd:
! - 
!>
  IMPLICIT NONE
  INTEGER, INTENT(IN)       ::  &
     &  jg                                       !< patch on which computation is performed
  TYPE(t_art_emiss_prescribed), INTENT(inout)  ::   &
     &  emiss_t                                  !< current prescribed emission dataset to be read
  TYPE(t_xml_file), INTENT(IN)  ::  &
     &  tixi_file_emiss                          !< XML file with emission dataset metadata
  !local variables
  LOGICAL                                  :: &
    &  l_exist                                   !< flag if file exists
  CHARACTER(LEN = max_datetime_str_len) ::  &
    &  datetime_str                              !< string of a datetime object
  CHARACTER(LEN=IART_PATH_LEN) ::  &
    &  emissiondataset
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo

  art_atmo => p_art_data(jg)%atmo


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Sanity check
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! search for the files corresponding to first and last date 
  ! and give an error message if non-existing
  CALL datetimeToPosixString(emiss_t%first_date,datetime_str,"%Y-%m-%d-%H")

  CALL art_create_filenames(jg,TRIM(emiss_t%path),emiss_t%iType_data, &
      &                       emissiondataset,.TRUE.,TRIM(datetime_str) )
  INQUIRE(FILE = TRIM(emissiondataset), EXIST=l_exist)
  IF (.NOT. l_exist) THEN
    CALL finish('mo_art_emiss_state:art_init_emission_dataset',  &
            &   'File of start date could not be found: '//TRIM(emissiondataset) &
            & //'. Check (1) link to "emissions" in'//TRIM(art_config(jg)%cart_input_folder) & 
            & //', (2) type of emission if it is one of ' &
            & //TRIM(art_get_all_emission_elements()) &
            & //' (also in the tracer xml file), (3) '  &
            & //'species, (4) inventory, (5) resolution ' &
            & //'and (6) start_date for dataset with id '//TRIM(emiss_t%id) &
            & //' in '//TRIM(tixi_file_emiss%path))
  END IF


  CALL datetimeToPosixString(emiss_t%last_date,datetime_str,"%Y-%m-%d-%H")

  CALL art_create_filenames(jg,TRIM(emiss_t%path),emiss_t%iType_data, &
      &                       emissiondataset,.TRUE.,TRIM(datetime_str) )

  INQUIRE(FILE = TRIM(emissiondataset), EXIST=l_exist)
  IF (.NOT. l_exist) THEN
    CALL finish('mo_art_emiss_state:art_init_emission_dataset',  &
            &   'File of end date could not be found: '//TRIM(emissiondataset) &
            & //'. Check end_date for dataset with id '//TRIM(emiss_t%id)//' in ' &
            & //TRIM(tixi_file_emiss%path))
  END IF

  emiss_t%first_date_is_leap_year = (getNoOfDaysInYearDateTime(emiss_t%first_date) == 366)
  emiss_t%last_date_is_leap_year  = (getNoOfDaysInYearDateTime(emiss_t%last_date) == 366)

  NULLIFY(art_atmo)
END SUBROUTINE art_init_emission_dataset
!!
!!--------------------------------------------------------------------------------------
!!
SUBROUTINE art_aerosol_assign_emission2tracer(p_mode_state, artconf, this_list, & 
  &                                           p_art_data, caero_emiss_xml_file)
!<
! SUBROUTINE art_aerosol_assign_emission2tracer
! This subroutine creates a link between emission routines and tracers.
! Part of Module: mo_art_emiss_state
! Author: Daniel Rieger, KIT
! Initial Release: 2017-04-26
! Modifications:
! 2018-10-25: Simon Gruber, KIT
! - modified for allowing flexible definiton of species, modes
!   and initial emission diameters
! yyyy-mm-dd:
! - 
!>
  TYPE(t_art_config), INTENT(IN)   :: &
    &  artconf                        !< ART configuration state
  TYPE(t_var_list_ptr),INTENT(IN)  :: &
    &  this_list                      !< current list: prognostic
  TYPE(t_art_data), TARGET, INTENT(inout) :: &
    &  p_art_data                     !< ART data container
  TYPE(t_mode_state),INTENT(inout) :: &
    &  p_mode_state
  CHARACTER(LEN=*),INTENT(IN)      :: &
    &  caero_emiss_xml_file           !< XML file containing information about emission routines
! Local variables
  TYPE(t_mode),POINTER             :: &
    &  current
  TYPE(t_mode_list),TARGET,SAVE    :: &
    &  newlist
  INTEGER                          :: &
    &  nemis                          !< Number of emission routines
  CHARACTER(LEN=1023)              :: &
    &  dummyStr                       !<
   
  NULLIFY(current)

  CALL new_mode_list(newlist)
  p_art_data%tracer2aeroemiss%e2t_list => newlist

  SELECT CASE (TRIM(caero_emiss_xml_file))
    CASE('')
      ! Do nothing, p_art_data%tracer2aeroemiss%lisinit remains .FALSE.
    CASE DEFAULT
      p_art_data%tracer2aeroemiss%lisinit = .TRUE.
      ! read emission routine specifications from xml
      CALL art_init_aerosol_emissions_xml(p_art_data, caero_emiss_xml_file, nemis)

      ! fill emission structure with metadata from tracers and modes
      CALL art_init_aerosol_assign_tracer(this_list, p_art_data, p_mode_state, nemis)

      ! calculate weights and print meta data
      current=>p_art_data%tracer2aeroemiss%e2t_list%p%first_mode
      DO WHILE(ASSOCIATED(current))
        SELECT TYPE(this=>current%fields)
          TYPE IS(t_art_emiss2tracer)
            SELECT CASE(this%name)
              CASE('seas_martensson','seas_monahan','seas_smith', &
                &  'seas_mode1', 'seas_mode2', 'seas_mode3')
                CALL this%calc_weights
              CASE DEFAULT
                this%weight(:,:) = 1.0_wp
            END SELECT

            CALL message(routine,TRIM(this%name))
            WRITE(dummyStr,*) this%nmodes, this%ntr
            CALL message(routine,'dims'//TRIM(dummyStr))
            WRITE(dummyStr,*) this%itr0
            CALL message(routine,'itr0'//TRIM(dummyStr))
            WRITE(dummyStr,*) this%itr3
            CALL message(routine,'itr3'//TRIM(dummyStr))
            WRITE(dummyStr,*) this%rho
            CALL message(routine,'rho'//TRIM(dummyStr))
            WRITE(dummyStr,*) this%sigmag
            CALL message(routine,'sigmag'//TRIM(dummyStr))
            WRITE(dummyStr,*) this%dg0
            CALL message(routine,'dg0'//TRIM(dummyStr))
            WRITE(dummyStr,*) this%dg3
            CALL message(routine,'dg3'//TRIM(dummyStr))
        END SELECT
        current=>current%next_mode
      END DO

      NULLIFY(current)
  END SELECT

END SUBROUTINE art_aerosol_assign_emission2tracer
!!
!!------------------------------------------------------------------------------------
!!
SUBROUTINE art_init_aerosol_assign_tracer(this_list, p_art_data, p_mode_state, nemis)
!<
! SUBROUTINE art_init_aerosol_assign_tracer
! This subroutine creates a link between emission routines and tracers.
! Part of Module: mo_art_emiss_state
! Author: Daniel Rieger, KIT
! Initial Release: 2017-04-26
! Modifications:
! 2018-10-25: Simon Gruber, KIT
! - modified for allowing flexible definiton of species, modes 
!   and initial emission diameters
! 2020-mm-dd: Sven Werchner, KIT
! - added additional flexibility and generic behaviour by using a linked-list
! - redesigned from "tracers are emitted by routine" to "routines emit tracers"
! - 
!>
  TYPE(t_var_list_ptr),INTENT(IN) :: &
    &  this_list                       !< current list: prognostic
  TYPE(t_art_data), TARGET, INTENT(inout) :: &
    &  p_art_data                      !< ART data container
  TYPE(t_mode_state),INTENT(inout) :: &
    &  p_mode_state
  INTEGER, INTENT(IN)            :: &
    &  nemis                           !< Number of emission routines
! Local variables
  TYPE(t_art_emiss2tracer),POINTER :: &
    &  this                            !< Current emiss2tracer dictionary item
  TYPE(t_var_metadata), POINTER   :: &
    &  info                            !< Returns reference to tracer metadata of current element
  TYPE(t_var_metadata_dynamic), POINTER :: &
    &  info_dyn                        !< Returns reference to dynamic tracer metadata of 
                                       !   current element
  TYPE(t_mode),POINTER             :: &
    &  this_mode,                     &  !< Pointer to current mode
    &  current_e2t
  INTEGER                         :: &
    &  idx,                          & !< Loop index
    &  imodes,                       & !< Loop index     
    &  itr,                          & !< Loop index  
    &  iemis,                        & !< Loop index
    &  itr_emis,                     & !< Loop index
    &  imode_emis,                   & !< Loop index
    &  isub,                         & !< Loop index
    &  iv,                           & !< Loop index
    &  loopCounter,                  & !< counts complete loop_runthroughs
    &  ierror                          !< Error return value
  CHARACTER(LEN=IART_VARNAMELEN)  :: &
    !&  caeroemiss,                   & !< Emission scheme associated to tracer
    &  substance,                    & !< Name of emitted substance
    &  mode_str,                     & !< Metadata mode
    &  new_modename,                 & !< temporary variable
    &  tracer_name                     !< Name of tracer
  REAL(wp)                        :: &
    &  log10dg0,log10diamode,        & !< needed for mode-determination
    &  delta_min

  NULLIFY(current_e2t)
  NULLIFY(info_dyn)
  NULLIFY(info)
  
  DO loopCounter=1,2
    mode_str = ''
    current_e2t => p_art_data%tracer2aeroemiss%e2t_list%p%first_mode
    DO WHILE (ASSOCIATED(current_e2t))
      SELECT TYPE(current_fields=>current_e2t%fields)
        TYPE IS(t_art_emiss2tracer)
          DO iemis=1,current_fields%nmodes
            log10dg0 = LOG10(current_fields%dg0(iemis,1))
            DO isub=1,current_fields%nsub
              SELECT CASE(loopCounter)
                CASE(1) ! preparation for modes and ntr determination
                  delta_min = 999.0_wp
                  new_modename = ''
                CASE(2) ! Allocating final structure if necessary
                  current_fields%ntr = MAXVAL(current_fields%ntrpermode)
                  IF(.NOT.ALLOCATED(current_fields%itr0))  &
                    & ALLOCATE(current_fields%itr0(current_fields%nmodes))
                  IF(.NOT.ALLOCATED(current_fields%itr3)) THEN
                    ALLOCATE(current_fields%itr3(current_fields%nmodes, &
                      &                           current_fields%ntr))
                    current_fields%itr3 = 0
                  END IF
                  IF(.NOT.ALLOCATED(current_fields%sigmag))                 &
                    & ALLOCATE(current_fields%sigmag(current_fields%nmodes, &
                    &                                 current_fields%ntr))
                  IF(.NOT.ALLOCATED(current_fields%molweight))                 & 
                    & ALLOCATE(current_fields%molweight(current_fields%nmodes, &
                    &                                   current_fields%ntr))
                  IF(.NOT.ALLOCATED(current_fields%weight))                 &
                    & ALLOCATE(current_fields%weight(current_fields%nmodes, &
                    &                                current_fields%ntr))
                CASE DEFAULT
                  ! Nothing to do
              END SELECT
              substance = current_fields%substance(isub)
              DO iv = 1, this_list%p%nvars

                info=>this_list%p%vl(iv)%p%info
                info_dyn=>this_list%p%vl(iv)%p%info_dyn
                IF (info_dyn%tracer%lis_tracer) THEN
                  SELECT TYPE(tracer_info => info_dyn%tracer)
                    CLASS IS(t_aero_meta)
                      tracer_name = TRIM(tracer_info%name)
                      IF(TRIM(tracer_info%substance)==TRIM(substance)) THEN
                        current_fields%lcalcemiss = .TRUE.
                        ! temporarily save name of corresponding mode
                        mode_str=TRIM(tracer_info%mode)                 
                        ! find modes' metadata 
                        NULLIFY(this_mode)
                        this_mode => p_mode_state%p_mode_list%p%first_mode
                        DO WHILE(ASSOCIATED(this_mode))
                          SELECT TYPE(this_fields => this_mode%fields)
                            CLASS IS(t_fields_2mom)
                              IF (TRIM(this_fields%name) == TRIM(mode_str)) THEN 
                                ! found mode containing current tracer
                                SELECT CASE(loopCounter)
                                  CASE(1) !determining possible modes and ntr
                                    log10diamode = LOG10(this_fields%info%diameter_ini_nmb)
                                    IF(ABS(log10diamode-log10dg0) < delta_min) THEN
                                      delta_min = ABS(log10diamode-log10dg0)
                                      new_modename = TRIM(mode_str)
                                    END IF
                                  CASE(2) !filling final structure
                                    IF(INDEX(TRIM(current_fields%modenames(iemis)),               &
                                      &                           TRIM(this_fields%name))>0) THEN
                                      DO itr=1,current_fields%ntr
                                        IF(current_fields%itr3(iemis,itr)==0) THEN
                                          current_fields%itr3(iemis,itr)=info%ncontained
                                          current_fields%itr0(iemis)=this_fields%itr0
                                          ! set data
                                          CALL info_dyn%tracer%opt_meta%get('mol_weight', &
                                            &         current_fields%molweight(iemis,itr), ierror)
                                          IF (ierror /= SUCCESS) CALL finish(TRIM(routine)// &
                                            &   ':art_aerosol_assign_emission2tracer',       &
                                            &   'Required metadata mol weight missing: '//   &
                                            &   TRIM(tracer_name))
                                          !DEBUG
                                          IF (current_fields%itr3(iemis,itr) /= 0 .AND. & 
                                            & current_fields%itr0(iemis) /= 0) THEN
                                            WRITE (message_text,*) 'ART: Emissions scheme '      // &
                                              &        TRIM(current_fields%name)//' assigned to '// &
                                              &        TRIM(tracer_name)
                                            CALL message(TRIM(routine)// &
                                              &    ':art_aerosol_assign_emission2tracer', message_text)
                                          ELSE
                                            CALL finish(TRIM(routine)//                         &
                                              &  ':art_aerosol_assign_emission2tracer',         &
                                              &  'Association of tracer '//TRIM(tracer_name)//  &
                                              &  ' to '//TRIM(current_fields%name)//' not successful.')
                                          ENDIF
                                          !END DEBUG
                                          EXIT
                                        ELSE IF(current_fields%itr3(iemis,itr) ==  &
                                          &     info%ncontained) THEN
                                          WRITE(message_text,*) 'ART: Emissions scheme already'// &
                                            &                   ' assigned.'
                                          CALL message(TRIM(routine)// &
                                            & ':art_aerosol_assign_emission2tracer', message_text)
                                          EXIT
                                        END IF
                                      END DO
                                    END IF
                                  CASE DEFAULT
                                  ! Nothing to do - should never be reached
                                END SELECT !loopCounter
                              ENDIF ! TRIM(this_fields%name) == TRIM(mode_str)
                            CLASS DEFAULT
                            ! Nothing to do
                          END SELECT !this_fields
                          this_mode => this_mode%next_mode
                        END DO !ASSOCIATED(this_mode)
                      END IF ! tracer_name==TRIM(substance)
                    CLASS DEFAULT
                    ! Nothing to do
                  END SELECT !tracer_info
                END IF !lis_tracer
              END DO !iv
              SELECT CASE(loopCounter)
                CASE(1)
                  IF(INDEX(current_fields%modenames(iemis),TRIM(new_modename))==0)&
                    & current_fields%modenames(iemis) = TRIM(new_modename)        &
                    &             //","//TRIM(current_fields%modenames(iemis))

                  current_fields%ntrpermode(iemis) =  &
                    &                 MAX(current_fields%ntrpermode(iemis)+1,1)
                CASE DEFAULT
                ! Nothing to do
              END SELECT !loopCounter
            END DO ! nsub
          END DO !nmodes
        CLASS DEFAULT
        ! Nothing to do
      END SELECT !current_fields
      current_e2t => current_e2t%next_mode
    END DO !ASSOCIATED(current_e2t)
  END DO !loopCounter
  
  NULLIFY(this_mode)
  NULLIFY(current_e2t)
  NULLIFY(info_dyn)
  NULLIFY(info)

END SUBROUTINE art_init_aerosol_assign_tracer
!!
!!------------------------------------------------------------------------------------
!!
SUBROUTINE art_init_aerosol_emissions_xml(p_art_data, caero_emiss_xml_file, nemis)
!<
! SUBROUTINE art_init_aerosol_emissions_xml
! This subroutine fills the structure for aerosol emission routines
! Part of Module: mo_art_emiss_state
! Author: Simon Gruber, KIT
! Initial Release: 2018-10-09
! Modifications:
! yyyy-mm-dd:
! -
!>
  TYPE(t_art_data), TARGET, INTENT(INOUT) :: &
    &  p_art_data                     !< ART data container
  CHARACTER(LEN=*),INTENT(IN)    :: &
    &  caero_emiss_xml_file           !< XML file containing point source information
  INTEGER, INTENT(INOUT)         :: &
    &  nemis                          !< Number of emission routines
! Local variables
  !TYPE(t_art_emiss2tracer), POINTER :: &
  TYPE(t_mode),POINTER :: &
    &  this                           !< Current emiss2tracer dictionary item
  TYPE(t_xml_file)               :: &
    &  aero_emiss_xmlfile             !< point source XML file
  TYPE(t_key_value_store)        :: &
    &  storage                        !< pointer to meta storage
  INTEGER                        :: &
    &  idx,                         & !< Loop index
    &  imodes,                      & !< Loop index
    &  nmodes,                      & !< Number of modes served by current emission routine
    &  ierror
  REAL(wp)                       :: &
    &  dg0,                         & !<
    &  dg3,                         & !<
    &  sigmag,                      & !<
    &  rho                            !<
  CHARACTER(LEN=IART_VARNAMELEN) :: &
    &  idx_str,                     & !<
    &  imodes_str,                  & !<
    &  emiss_name                     !< Name of source
  CHARACTER(:), ALLOCATABLE      :: &
    &  c_tmp
  ! temporary - maybe there is a better solution
  INTEGER           :: &
    &  idims_dummy(3)
    
  idims_dummy(:) = 1

  CALL art_open_xml_file(TRIM(caero_emiss_xml_file), aero_emiss_xmlfile)

  ! read xml file
  CALL art_get_childnumber_xml(aero_emiss_xmlfile, '/emiss/', nemis)
  WRITE (message_text,'(A31,I4)') 'ART: Number of aerosol emission routines: ',nemis
  CALL message (TRIM(routine)//':art_init_aero_emiss', message_text)

  DO idx = 1, nemis
    CALL storage%init(.FALSE.)

    ! read name and metadata of current emission routine in XML file
    WRITE(idx_str,'(I3)') idx
    CALL art_read_elements_xml(aero_emiss_xmlfile,'/emiss/*['//TRIM(ADJUSTL(idx_str))//']/', &
      &                        storage)
    CALL storage%get('name',c_tmp,ierror)
    IF(ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_init_aero_emiss',&
                            &         'name not found in XML file.')
    WRITE(emiss_name,'(A)') c_tmp
    CALL storage%get('nmodes',nmodes,ierror)
    IF(ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_init_aero_emiss',&
                            &         'nmodes not found in XML file.')
    ! rho is required in every case
    CALL storage%get('rho',rho,ierror)
    IF(ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_init_aero_emiss',&
                            &         'rho not found in XML file.')

    ! Prepare and fill corresponding emission structure
    NULLIFY(this)
    CALL append_mode(p_art_data%tracer2aeroemiss%e2t_list,this,IART_EMISS2TRACER)
    !fields allocation done here, since it was not possible in append_mode->create_mode
    ALLOCATE(t_art_emiss2tracer   :: this%fields)
    CALL this%fields%create(TRIM(emiss_name),idims_dummy) !emiss_name='seas','dust','volc',...

    IF (ASSOCIATED(this)) THEN
      CALL this%fields%set_meta(storage)
    ENDIF ! ASSOCIATED(this)
    NULLIFY(this)
    CALL storage%destruct
  ENDDO ! idx
  CALL art_close_xml_file(aero_emiss_xmlfile)

END SUBROUTINE art_init_aerosol_emissions_xml
!!
!!--------------------------------------------------------------------------------------
!!
END MODULE mo_art_emiss_state
