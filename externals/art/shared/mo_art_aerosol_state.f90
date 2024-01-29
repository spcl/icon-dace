!
! mo_art_aerosol_state
! This module provides initialization for the aerosol structures
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

MODULE mo_art_aerosol_state
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_exception,                     ONLY: message, finish, message_text
  USE mo_var_list,                      ONLY: t_var_list_ptr
  USE mo_var_metadata_types,            ONLY: t_var_metadata, t_var_metadata_dynamic
  USE mo_tracer_metadata_types,         ONLY: t_aero_meta
  USE mo_key_value_store,               ONLY: t_key_value_store
  USE mo_impl_constants,                ONLY: SUCCESS
  USE mo_math_constants,                ONLY: pi
  USE mo_var_list,                      ONLY: get_tracer_info_dyn_by_idx
! ART
  USE mo_art_modes_linked_list,         ONLY: t_mode_state, mode_lists,        &
                                          &   n_mode_lists, max_mode_lists,    &
                                          &   new_mode_list, t_mode, append_mode
  USE mo_art_modes,                     ONLY: t_fields_1mom, t_fields_2mom,    &
                                          &   t_fields_radio, t_fields_volc
  USE mo_art_read_xml,                  ONLY: t_xml_file, art_open_xml_file,   &
                                          &   art_close_xml_file,              &
                                          &   art_get_childnumber_xml,         &
                                          &   art_read_elements_xml
  USE mo_art_impl_constants,            ONLY: IART_VARNAMELEN,                 &
                                          &   IART_MODE_RADIO, IART_MODE_VOLC, &
                                          &   IART_MODE_POLL, IART_MODE_2MOM
  USE mo_art_string_tools,              ONLY: key_value_storage_as_string

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_aerosol_state'

  PUBLIC :: art_aerosol_state
  PUBLIC :: art_assign_tracers_to_modes
  PUBLIC :: set_meta_init_nmb_mass_conc

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_aerosol_state(p_mode_state, list_suffix, xml_file, coag_xml, idims)
!<
! SUBROUTINE art_aerosol_state
! This subroutine provides the initialization routines for aerosol structures.
! Part of Module: mo_art_aerosol_state
! Author: Daniel Rieger, KIT
! Initial Release: 2016-08-19
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  TYPE(t_mode_state),INTENT(inout) :: &
    &  p_mode_state
  CHARACTER(LEN=*),INTENT(in)      :: &
    &  list_suffix,                   & !< Suffix of the modes list name. (=jg for ICON-ART)
    &  xml_file,                      & !< File name and path to XML file for mode initialization
    &  coag_xml                         !< File name and path to XML file for coagulation initialization
  INTEGER,INTENT(in)               :: &
    &  idims(3)                         !< Dimensions to allocate fields
!Local variables
  TYPE(t_xml_file)                 :: &
    &  tixi_file                        !< Filename and TIXI handle of XML file
  TYPE(t_key_value_store)          :: &
    &  key_value_store                  !< Metadata for mode
  TYPE(t_mode),POINTER             :: &
    &  this_mode,                     & !< Pointer to current mode
    &  partner_mode,                  & !< Pointer to partner mode (coagulation)
    &  target_mode                      !< Pointer to target mode (coagulation)
  INTEGER                          :: &
    &  i, idx_mode_xml,               & !< Counter
    &  nmodes_xml,                    & !< Number of modes in XML file
    &  nmodes_add,                    & !< number of partnermodes in XML file
    &  IART_MOMENT,                   & !< Mode type (IART_MODE_VOLC, IART_MODE_2MOM, ...)
    &  ierror,                        & !< Error variable
    &  thisIdx, partnerIdx              !< counter for coagulation-loop
  CHARACTER(LEN=IART_VARNAMELEN)   :: &
    &  mode_list_name                   !< name of the mode list to be created
  CHARACTER(:), ALLOCATABLE        :: &
    &  modename,                      & !< Name of mode
    &  modekind,                      & !< Kind of mode (2mom, 1mom, ...)
    &  c_targetmode
  CHARACTER(LEN=3)                 :: &
    &  idx_mode_str                     !< idx_mode_xml as character
  TYPE(t_fields_2mom)              :: &
    &  this_fields                      !< fields of this_mode
  !TYPE(t_fields_2mom),POINTER      :: &
  !  &  partner_fields,                & !< fields of partner_mode (coagulation)
  !  &  target_fields                    !< fields of target_mode (coagulation)
  LOGICAL                          :: &
    &  lfound                           !< Flag for searching modes
  TYPE(t_mode),POINTER             :: &
    &  that_mode                        !< Pointer for mode shifting

  ! ----------------------------------
  ! --- Create mode structure
  ! ----------------------------------

  mode_list_name = 'ART_MODE_'//TRIM(list_suffix)

  WRITE (message_text,*) 'ART: Creating new mode list '//TRIM(mode_list_name)
  CALL message (TRIM(routine)//':art_aerosol_state', message_text)

  ! look, if name exists already in list
  DO i = 1, n_mode_lists
    IF (TRIM(mode_lists(i)%p%name) == TRIM(mode_list_name)) THEN
      CALL finish(TRIM(routine)//':art_aerosol_state',          &
        &         'mode_list '//TRIM(mode_list_name)//' already used.')
    ENDIF
  ENDDO

  IF(.NOT. ASSOCIATED(p_mode_state%p_mode_list)) THEN
    n_mode_lists = n_mode_lists + 1
    IF (n_mode_lists > max_mode_lists) THEN
      CALL finish(TRIM(routine)//':art_aerosol_state',          &
        &         'mode_lists container overflow, increase max_mode_lists.')
    ENDIF
  ELSE
    CALL finish(TRIM(routine)//':art_aerosol_state',          &
      &         'p_mode_state%p_mode_list already associated.')
  ENDIF

  ! Eventually create the mode list
  CALL new_mode_list(mode_lists(n_mode_lists))

  ! connect anchor and backbone by referencing and give it a name
  p_mode_state%p_mode_list => mode_lists(n_mode_lists)
  p_mode_state%p_mode_list%p%name = TRIM(mode_list_name)

  ! ----------------------------------
  ! --- Initialize mode structure with XML file
  ! ----------------------------------

  CALL art_open_xml_file(TRIM(xml_file),tixi_file)

  CALL art_get_childnumber_xml(tixi_file, '/modes/', nmodes_xml)
  WRITE (message_text,*) 'ART: Appending modes to '//TRIM(mode_list_name)
  CALL message (TRIM(routine)//':art_aerosol_state', message_text)

  DO idx_mode_xml = 1, nmodes_xml

    WRITE(idx_mode_str,'(I3)') idx_mode_xml
    CALL art_read_elements_xml(tixi_file,'/modes/*['//TRIM(ADJUSTL(idx_mode_str))//']/',key_value_store)
    CALL key_value_storage_as_string(key_value_store,'name',modename,ierror)

    IF(ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_aerosol_state',    &
                            &         'Mode name not available in meta data.')
    CALL key_value_storage_as_string(key_value_store,'kind',modekind,ierror)

    IF(ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_aerosol_state',    &
                            &         'Mode kind not available in meta data.')

    IART_MOMENT = 0
    IF(TRIM(modekind)=='2mom')  IART_MOMENT = IART_MODE_2MOM
    IF(TRIM(modekind)=='radio') IART_MOMENT = IART_MODE_RADIO
    IF(TRIM(modekind)=='volc')  IART_MOMENT = IART_MODE_VOLC
    IF(TRIM(modekind)=='poll')  IART_MOMENT = IART_MODE_POLL
    IF (IART_MOMENT == 0) THEN
      CALL finish(TRIM(routine)//':art_aerosol_state', &
        &         'Mode type '//TRIM(modekind)//' unknown.')
    ENDIF

    WRITE (message_text,*) 'ART: Appending mode '//TRIM(modename)//'of kind '//TRIM(modekind)//'.'
    CALL message (TRIM(routine)//':art_aerosol_state', message_text)

    CALL append_mode(p_mode_state%p_mode_list, this_mode, IART_MOMENT)

    ! Create structure (set name, allocate, set initial undefined values)
    CALL this_mode%fields%create(TRIM(modename), idims)

    ! Write metadata from key_value_store to this_mode%fields
    CALL this_mode%fields%set_meta(key_value_store)         

    ! Memory allocated for this_mode is accessible via p_mode_list,
    ! so we can nullify the pointer in order to reuse it
    NULLIFY(this_mode)

    CALL key_value_store%destruct
  ENDDO !idx_mode_xml

  CALL art_close_xml_file(tixi_file)
  
  
  ! ----------------------------------
  ! --- Initialize mode coagulation: connect modes ("this_mode + partner_mode = target_mode")
  ! ----------------------------------
  ! - NOTE: '=' instead of '=>' at field-assignments needs to be checked for functionality
  ! ----------------------------------
  ! coagulation-xmlfile
  IF(TRIM(coag_xml)/='') CALL art_open_xml_file(TRIM(coag_xml),tixi_file)
  this_mode => p_mode_state%p_mode_list%p%first_mode
  thisIdx = 1
  DO WHILE(ASSOCIATED(this_mode))
    SELECT TYPE(this_fields => this_mode%fields)
      CLASS IS(t_fields_2mom)
        ! Is coagulation enabled for this mode?
        IF(this_fields%do_coag <= 0) THEN
          this_mode => this_mode%next_mode
          CYCLE
        END IF
        IF(TRIM(coag_xml)=='') THEN
          CALL finish (TRIM(routine)//':art_aerosol_state', 'coagulation tag found, but no '  &
            &                                               //'coagulation-XML-File given')
        END IF
        ! Search for partner modes
        ! Get XML-Data for coagulation
        CALL art_read_elements_xml(tixi_file,'/coagulate/smallmode[@id="' &
          &                        //TRIM(ADJUSTL(this_fields%name))//'"]/', key_value_store)
        CALL key_value_store%get('nmodes',nmodes_add, ierror)
        IF(ierror/=0) THEN
          CALL finish (TRIM(routine)//':art_aerosol_state', 'coagulation enabled without specification in XML-File')
        END IF
        this_fields%coag_util%n_modes = this_fields%coag_util%n_modes + nmodes_add
        ! Allocating
        IF (.NOT. ALLOCATED(this_fields%coag_util%p_coag)) THEN
          ALLOCATE(this_fields%coag_util%p_coag(this_fields%coag_util%n_modes))
        END IF
        ! Setting partner and target mode
        partner_mode => p_mode_state%p_mode_list%p%first_mode
        partnerIdx = 1
        DO WHILE(ASSOCIATED(partner_mode))
          SELECT TYPE(partner_fields => partner_mode%fields)
            CLASS IS(t_fields_2mom)
              ! Is coagulation enabled for this mode?
              IF(partner_fields%do_coag <= 0) THEN
                partner_mode => partner_mode%next_mode
                CYCLE
              END IF
              ! Check if already dealt with
              IF ( partnerIdx < thisIdx ) THEN
                this_fields%coag_util%p_coag(partnerIdx)%p_coagulateWith => partner_mode%fields
                this_fields%coag_util%p_coag(partnerIdx)%p_coagulateTo   => partner_fields%coag_util%p_coag(thisIdx)%p_coagulateTo
                ALLOCATE(this_fields%coag_util%p_coag(partnerIdx)%coagcoeff0(idims(1),idims(2),idims(3)))
                this_fields%coag_util%p_coag(partnerIdx)%coagcoeff0(:,:,:) = 0.0_wp
                ALLOCATE(this_fields%coag_util%p_coag(partnerIdx)%coagcoeff3(idims(1),idims(2),idims(3)))
                this_fields%coag_util%p_coag(partnerIdx)%coagcoeff3(:,:,:) = 0.0_wp
                WRITE(message_text,*) 'ART: Mode '//TRIM(this_fields%name)//                      &
                &                   ' will be able to coagulate with mode '//                     &
                &                   TRIM(this_fields%coag_util%p_coag(partnerIdx)%p_coagulateWith%name)//&
                &                   ' to mode '//                                                 &
                &                   TRIM(this_fields%coag_util%p_coag(partnerIdx)%p_coagulateTo%name)//'.'
                CALL message (TRIM(routine)//':art_aerosol_state', message_text)
                partnerIdx = partnerIdx + 1
                partner_mode => partner_mode%next_mode
                CYCLE
              END IF
              CALL key_value_storage_as_string(key_value_store,'bigmode_'  &
                &                       //TRIM(ADJUSTL(partner_fields%name)),c_targetmode,ierror)
              IF(ierror/=0) THEN
                partner_mode => partner_mode%next_mode
                CYCLE
              ENDIF
              this_fields%coag_util%p_coag(partnerIdx)%p_coagulateWith => partner_mode%fields
              ALLOCATE(this_fields%coag_util%p_coag(partnerIdx)%coagcoeff0(idims(1),idims(2),idims(3)))
              this_fields%coag_util%p_coag(partnerIdx)%coagcoeff0(:,:,:) = 0.0_wp
              ALLOCATE(this_fields%coag_util%p_coag(partnerIdx)%coagcoeff3(idims(1),idims(2),idims(3)))
              this_fields%coag_util%p_coag(partnerIdx)%coagcoeff3(:,:,:) = 0.0_wp
              IF ( partnerIdx > thisIdx ) THEN
                partner_fields%coag_util%n_modes = partner_fields%coag_util%n_modes + 1
              END IF
              target_mode => p_mode_state%p_mode_list%p%first_mode
              DO WHILE(ASSOCIATED(target_mode))
                SELECT TYPE(target_fields => target_mode%fields)
                  CLASS IS(t_fields_2mom)
                    IF(TRIM(target_fields%name) == TRIM(c_targetmode)) THEN
                      ! target mode found, set pointer to it
                      this_fields%coag_util%p_coag(partnerIdx)%p_coagulateTo => target_mode%fields
                      EXIT
                    ENDIF
                  CLASS DEFAULT
                    ! Nothing to do
                END SELECT
                target_mode => target_mode%next_mode
              ENDDO
              WRITE(message_text,*) 'ART: Mode '//TRIM(this_fields%name)//                        &
                &                   ' will be able to coagulate with mode '//                     &
                &                   TRIM(this_fields%coag_util%p_coag(partnerIdx)%p_coagulateWith%name)//&
                &                   ' to mode '//                                                 &
                &                   TRIM(this_fields%coag_util%p_coag(partnerIdx)%p_coagulateTo%name)//'.'
              CALL message (TRIM(routine)//':art_aerosol_state', message_text)
              partnerIdx = partnerIdx + 1
            CLASS DEFAULT
              ! Nothing to do
          END SELECT
          partner_mode => partner_mode%next_mode
        ENDDO !ASSOCIATED(partner_mode)
        thisIdx = thisIdx + 1
        ! <vv - needs to be adjusted>
        !IF (this_fields%do_coag .AND. &
        !  & .NOT. ASSOCIATED(this_fields%coag_util%p_coagulateWith(i))) THEN
        !  CALL finish(TRIM(routine)//':art_aerosol_state', &
        !    &           'No partner modes found for '//TRIM(this_fields%name)//'.')
        !ENDIF
      CLASS DEFAULT
        ! Nothing to do
    END SELECT
 
    this_mode => this_mode%next_mode
  ENDDO 
  IF(TRIM(coag_xml)/='') CALL art_close_xml_file(tixi_file)

  ! ----------------------------------
  ! --- Initialize mode shifting: connect modes
  ! ----------------------------------
  this_mode => p_mode_state%p_mode_list%p%first_mode
  DO WHILE(ASSOCIATED(this_mode))
    SELECT TYPE(this_fields => this_mode%fields) 
      CLASS IS(t_fields_2mom)
        ! Search for mode to shift to
        that_mode => p_mode_state%p_mode_list%p%first_mode
        DO WHILE(ASSOCIATED(that_mode))
          ! Mode shift to larger mode
          IF(this_fields%shift2larger%l_do_shift) THEN
            IF (TRIM(that_mode%fields%name) == TRIM(this_fields%shift2larger%shift2name)) THEN
              ! Mode to shift found, set pointer to it
              IF(ASSOCIATED(this_fields%p_shift2larger)) CALL finish(TRIM(routine)//':art_aerosol_state', &
                                                           &         'Second mode to shift larger found.')
              this_fields%p_shift2larger => that_mode%fields
              WRITE (message_text,*) 'ART: Will shift mode '//TRIM(this_fields%name)// &
                &                    ' to larger mode '//TRIM(this_fields%p_shift2larger%name)//'.'
              CALL message (TRIM(routine)//':art_aerosol_state', message_text)
            ENDIF
          ENDIF
          ! Mode shift to mixed mode
          IF(this_fields%shift2mixed%l_do_shift) THEN
            IF (TRIM(that_mode%fields%name) == TRIM(this_fields%shift2mixed%shift2name)) THEN
              ! Mode to shift found, set pointer to it
              IF(ASSOCIATED(this_fields%p_shift2mixed)) CALL finish(TRIM(routine)//':art_aerosol_state', &
                                                          &         'Second mode to shift mixed found.')
              this_fields%p_shift2mixed => that_mode%fields
              WRITE (message_text,*) 'ART: Will shift mode '//TRIM(this_fields%name)//&
                 &                    ' to mixed mode '//TRIM(this_fields%p_shift2mixed%name)//'.'
              CALL message (TRIM(routine)//':art_aerosol_state', message_text)
            ENDIF
          ENDIF
          that_mode => that_mode%next_mode
        ENDDO

        IF (this_fields%shift2larger%l_do_shift .AND. &
          & .NOT. ASSOCIATED(this_fields%p_shift2larger)) THEN
          CALL finish(TRIM(routine)//':art_aerosol_state', &
            &           'No larger mode to shift found for '//TRIM(this_fields%name)//'.')
        ENDIF
        IF (this_fields%shift2mixed%l_do_shift .AND. &
          & .NOT. ASSOCIATED(this_fields%p_shift2mixed)) THEN
          CALL finish(TRIM(routine)//':art_aerosol_state', &
            &           'No mixed mode to shift found for '//TRIM(this_fields%name)//'.')
        ENDIF
      CLASS DEFAULT
        ! Nothing to do
    END SELECT

    this_mode => this_mode%next_mode
  ENDDO

END SUBROUTINE art_aerosol_state
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_assign_tracers_to_modes(p_mode_state,this_list)
!<
! SUBROUTINE art_assign_tracers_to_modes
! This subroutine searches for tracers contained in the modes and sets according pointers
! to the mass mixing ratio fields and specific number fields
! Author: Daniel Rieger, KIT
! Initial Release: 2016-09-12
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  TYPE(t_mode_state),INTENT(inout) :: &
    &  p_mode_state                     !< mode container
  TYPE(t_var_list_ptr),INTENT(in)      :: &
    &  this_list                        !< current list: prognostic
! Local variables
  TYPE(t_mode),POINTER                  :: &
    &  this_mode                             !< Reference to current mode in list
  TYPE(t_fields_2mom),POINTER           :: &
    &  that_fields                           !< Reference to mode to shift to
  TYPE(t_var_metadata), POINTER         :: &
    &  info,                               & !< returns reference to dynamic tracer metadata of current element
    &  that_info                             !< returns reference to dynamic tracer metadata of current element
  TYPE(t_var_metadata_dynamic), POINTER :: &
    &  info_dyn,                           & !< returns reference to dynamic tracer metadata of current element
    &  that_info_dyn                         !< returns reference to dynamic tracer metadata of current element
  INTEGER                               :: &
    &  ierror,                             & !< Error return value
    &  inucleation,                        & !< Nucleation scheme
    &  ntracer_in_mode,                    & !< Number of tracer in current mode
    &  ntracer,                            & !< Number of tracer
    &  iv                                    !< loop index
  CHARACTER(LEN=2)                      :: &
    &  cntracer                              !< Number of tracer as character
  CHARACTER(LEN=50)                     :: &
    &  this_tracer_name                      !< short name of tracer

  this_mode => p_mode_state%p_mode_list%p%first_mode

  DO WHILE(ASSOCIATED(this_mode))
    ! First step: Find out how many tracer are contained in the mode
    ntracer_in_mode = 0
    ntracer         = 0
    DO iv = 1, this_list%p%nvars
      info_dyn=>this_list%p%vl(iv)%p%info_dyn
      IF (info_dyn%tracer%lis_tracer) THEN
        SELECT TYPE(tracer_info => info_dyn%tracer)
          CLASS IS(t_aero_meta)
            IF (TRIM(this_mode%fields%name) == TRIM(tracer_info%mode)) THEN
              ntracer_in_mode = ntracer_in_mode + 1
            ENDIF
          CLASS DEFAULT
            ! Nothing to do for other tracer
        END SELECT
      ENDIF
    ENDDO
    
    ! Write number of tracer to character for output statement
    IF (ntracer_in_mode .le. 9) THEN
      write(cntracer,"(I1)") ntracer_in_mode
      cntracer = '0'//cntracer
    ELSE
      write(cntracer,"(I2)") ntracer_in_mode
    ENDIF
    
    WRITE (message_text,*) 'ART: Mode '//TRIM(this_mode%fields%name)//' contains '//cntracer//' tracer.'
    CALL message (TRIM(routine)//':art_assign_tracers_to_modes', message_text)

    SELECT TYPE(fields => this_mode%fields)
!--------------------------------------------------------------------------------------------------
      CLASS IS(t_fields_2mom)
        ! Allocate mass mixing ratio pointer array. -1 as one of the tracer is 
        ! assumed to be the specific number of the mode
        ALLOCATE(fields%itr3(ntracer_in_mode-1))
        ALLOCATE(fields%linsol(ntracer_in_mode-1))
        ALLOCATE(fields%shift2mixed%itr3(2,ntracer_in_mode-1))
        ALLOCATE(fields%shift2mixed%itr0(2))
        fields%shift2mixed%njsp  = ntracer_in_mode-1
        ALLOCATE(fields%shift2larger%itr3(2,ntracer_in_mode-1))
        ALLOCATE(fields%shift2larger%itr0(2))
        fields%shift2larger%njsp = ntracer_in_mode-1
        ALLOCATE(fields%rho(ntracer_in_mode-1))
        ALLOCATE(fields%info%init_mass_conc(ntracer_in_mode-1))
        fields%ntr = ntracer_in_mode
        ! Connect
        DO iv = 1, this_list%p%nvars
          info    =>this_list%p%vl(iv)%p%info
          info_dyn=>this_list%p%vl(iv)%p%info_dyn
          IF (info_dyn%tracer%lis_tracer) THEN
            SELECT TYPE(tracer_info => info_dyn%tracer)
              CLASS IS(t_aero_meta)
                IF (TRIM(fields%name) == TRIM(tracer_info%mode)) THEN
                  IF (tracer_info%moment == 0) THEN
                    fields%itr0 = info%ncontained
                    WRITE (message_text,*) 'ART: '//TRIM(info%name)//' is number of '//TRIM(fields%name)
                    CALL message (TRIM(routine)//':art_assign_tracers_to_modes', message_text)
                  ENDIF
                  IF (tracer_info%moment == 3) THEN
                    ntracer = ntracer + 1
                    IF (ntracer > ntracer_in_mode-1) THEN
                      CALL finish(TRIM(routine)//':art_assign_tracers_to_modes', &
                        &         'Trying to add too much tracer to mode '//TRIM(fields%name))
                    ENDIF
                    fields%itr3(ntracer) =  info%ncontained
                    WRITE (message_text,*) 'ART: '//TRIM(info%name)//' is mass of '//TRIM(fields%name)
                    CALL message (TRIM(routine)//':art_assign_tracers_to_modes', message_text)
                    CALL info_dyn%tracer%opt_meta%get('rho',fields%rho(ntracer),ierror)
                    IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_assign_tracers_to_modes', &
                      &                                'Required metadata rho not found.')
                    fields%linsol(ntracer) = tracer_info%linsol
                     ! Nucleation initialization
                    CALL info_dyn%tracer%opt_meta%get('inucleation',inucleation,ierror)
                    IF (ierror == SUCCESS) THEN
                      WRITE (message_text,*) 'ART: Initializing nucleation for '//TRIM(info%name)//'.'
                      CALL message (TRIM(routine)//':art_assign_tracers_to_modes', message_text)
                      fields%info%inucleation = inucleation
                      fields%info%itrnucl     = info%ncontained
                    ENDIF
                    ! Condensation initialization
                    IF (fields%info%icondensation > 0) THEN
                      ! search for a tracer with so4 in name
                      IF (index(TRIM(info%name),'so4') > 0) THEN
                        fields%info%itrcond   = info%ncontained
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
              CLASS DEFAULT
                ! Nothing to do for other tracer
            END SELECT
          ENDIF
        ENDDO
        IF (ntracer /= ntracer_in_mode-1) THEN
          CALL finish(TRIM(routine)//':art_assign_tracers_to_modes', &
            &         'ntracer and ntracer_in_mode-1 do not match.')
        ENDIF
!--------------------------------------------------------------------------------------------------
      CLASS IS(t_fields_1mom)
        IF (ntracer_in_mode > 1) THEN
          CALL finish(TRIM(routine)//':art_assign_tracers_to_modes', &
            &         '1MOM-Mode '//TRIM(fields%name)//' contains more than one tracer.')
        ENDIF
        ! Connect
        DO iv = 1, this_list%p%nvars
          info    =>this_list%p%vl(iv)%p%info
          info_dyn=>this_list%p%vl(iv)%p%info_dyn
          IF (info_dyn%tracer%lis_tracer) THEN
            SELECT TYPE(tracer_info => info_dyn%tracer)
              CLASS IS(t_aero_meta)
                IF (TRIM(fields%name) == TRIM(tracer_info%mode)) THEN
                  fields%itr = info%ncontained
                  WRITE (message_text,*) 'ART: '//TRIM(info%name)//' is part of '//TRIM(fields%name)
                  CALL message (TRIM(routine)//':art_assign_tracers_to_modes', message_text)
                ENDIF
              CLASS DEFAULT
                ! Nothing to do for other tracer
            END SELECT
          ENDIF
        ENDDO
!--------------------------------------------------------------------------------------------------
      CLASS DEFAULT
        CALL finish(TRIM(routine)//':art_assign_tracers_to_modes', &
          &         'Unknown mode type for mode')
    END SELECT

    CALL this_mode%fields%print
    
    this_mode => this_mode%next_mode
  ENDDO


  this_mode => p_mode_state%p_mode_list%p%first_mode
  DO WHILE(ASSOCIATED(this_mode))
    SELECT TYPE(this_fields => this_mode%fields)
      CLASS IS(t_fields_2mom)
        IF (this_fields%shift2larger%l_do_shift) THEN
          SELECT TYPE(that_fields => this_fields%p_shift2larger)
            CLASS IS(t_fields_2mom)
              CALL connect_shifting(this_fields, that_fields, this_list, &
                &                   this_fields%shift2larger%itr3, this_fields%shift2larger%itr0)
            CLASS DEFAULT
              CALL finish(TRIM(routine)//':art_assign_tracers_to_modes',          &
                &         'p_shift2larger is not of type t_fields_2mom for'//TRIM(this_fields%name)//'.')
          END SELECT
        ENDIF
        IF (this_fields%shift2mixed%l_do_shift) THEN
          SELECT TYPE(that_fields => this_fields%p_shift2mixed)
            CLASS IS(t_fields_2mom)
              CALL connect_shifting(this_fields, that_fields, this_list, &
                &                   this_fields%shift2mixed%itr3, this_fields%shift2mixed%itr0)
            CLASS DEFAULT
              CALL finish(TRIM(routine)//':art_assign_tracers_to_modes',          &
                &         'p_shift2mixed is not of type t_fields_2mom for'//TRIM(this_fields%name)//'.')
          END SELECT
        ENDIF
      CLASS DEFAULT
        ! Nothing to do
    END SELECT
    this_mode => this_mode%next_mode
  ENDDO ! this_mode

END SUBROUTINE art_assign_tracers_to_modes
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE connect_shifting(shift_from, shift_to, this_list, itr3, itr0)
!<
! SUBROUTINE art_assign_tracers_to_modes
! This subroutine searches for tracers contained in the modes and sets according pointers
! to the mass mixing ratio fields and specific number fields
! Author: Daniel Rieger, KIT
! Initial Release: 2016-09-12
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  TYPE(t_fields_2mom), INTENT(inout) :: &
    &  shift_from,                  & !< Mode to shift from
    &  shift_to                       !< Mode to shift to
  TYPE(t_var_list_ptr),INTENT(in):: &
    &  this_list                      !< current list: prognostic
  INTEGER, INTENT(out)           :: &
    &  itr3(:,:), itr0(:)             !< Connectivity integer arrays
! Local variables
  TYPE(t_var_metadata_dynamic)   :: &
    &  info_dyn                       !< dynamic tracer metadata
  INTEGER                        :: &
    &  js1, js2,                    & !< Loop counter
    &  idx_us1, idx_us2               !< Index of first underscrore in string
  CHARACTER(LEN=IART_VARNAMELEN) :: &
    &  name1, name2,                & !< Complete tracer name
    &  subname1, subname2             !< Substring of tracer name 1 and 2

  itr0(1) = shift_from%itr0
  itr0(2) = shift_to%itr0

  DO js1 = 1, shift_from%ntr-1 ! Loop over mass species contained in mode
    CALL get_tracer_info_dyn_by_idx(this_list, shift_from%itr3(js1),info_dyn)
    !idx_us1 = INDEX(TRIM(info_dyn%tracer%name),'_')
    !subname1 = info_dyn%tracer%name(1:idx_us1-1)
    !name1    = TRIM(info_dyn%tracer%name)
    SELECT TYPE(tracer=>info_dyn%tracer)
      TYPE IS(t_aero_meta)
        subname1 = tracer%substance
      CLASS DEFAULT
        CALL finish(TRIM(routine)//':connect_shifting','shifting only applicable for t_aero_meta')
    END SELECT
    DO js2 = 1, shift_to%ntr-1 ! Loop over mass species contained in mode
      CALL get_tracer_info_dyn_by_idx(this_list, shift_to%itr3(js2),info_dyn)
      !idx_us2 = INDEX(TRIM(info_dyn%tracer%name),'_')
      !subname2 = info_dyn%tracer%name(1:idx_us2-1)
      !name2    = TRIM(info_dyn%tracer%name)
      SELECT TYPE(tracer=>info_dyn%tracer)
        TYPE IS(t_aero_meta)
          subname2 = tracer%substance
        CLASS DEFAULT
          CALL finish(TRIM(routine)//':connect_shifting','shifting only applicable for t_aero_meta')
      END SELECT
      IF(TRIM(subname1) == TRIM(subname2)) THEN
        itr3(1,js1) = shift_from%itr3(js1)
        itr3(2,js1) = shift_to%itr3(js2)
      ENDIF
    ENDDO
  ENDDO

! Checks
  DO js1 = 1, shift_from%ntr-1
    IF(.NOT.ANY(shift_from%itr3(js1) == itr3(1,:))) CALL finish(TRIM(routine)//':connect_shifting', &
      &                                                    'Connection not successful.')
    IF(.NOT.ANY(itr3(2,js1) == shift_to%itr3(:))) CALL finish(TRIM(routine)//':connect_shifting', &
      &                                                  'Connection not successful.')
  ENDDO

END SUBROUTINE connect_shifting
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE set_meta_init_nmb_mass_conc(p_mode_state)
!<
! SUBROUTINE set_meta_init_nmb_mass_conc
! This subroutine sets meta data for initial number and mass concentration in
! the
! fields_2mom meta data container. A fix number concentration of 100 #/kg per
! mode
! is prescribed. The corresponding mass mixing ratio is determined based on the
! amount of tracers per mode and the solubility class of the tracer.
! This meta data needs to be stored seperately, as the tracer initialization
! routine is not called during restart runs.
! Author: Lukas Muser, KIT
! Initial Release: 2021-02-26
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  TYPE(t_mode_state),INTENT(inout) :: &
    &  p_mode_state
  ! local variables
  TYPE(t_mode),POINTER             :: &
    &  current_mode
  INTEGER                          :: &
    &  ntr_insol, jsp
  REAL(wp)                         :: &
    &  v_conc
  REAL(wp), PARAMETER              :: &
    &  init_nmb = 100.0_wp

  ! loop through modes
  current_mode => p_mode_state%p_mode_list%p%first_mode
  DO WHILE(ASSOCIATED(current_mode))
    ! Select type of mode
    SELECT TYPE (fields => current_mode%fields)
      CLASS IS (t_fields_2mom)
        ! write init nmb conc to mode meta data
        fields%info%init_nmb_conc = init_nmb

        v_conc = pi / 6.0_wp * init_nmb / (fields%ntr-1) * fields%info%diameter_ini_nmb**3   &
          &    * fields%info%exp_aero**36 * 1.0E9_wp
        DO jsp=1,fields%ntr-1
          ! write init mass mixing ratio to mode meta data
          fields%info%init_mass_conc(jsp) = v_conc * fields%rho(jsp)
        ENDDO !jsp

        ! handel insoluble mode -> soluble fraction is initialized with 0
        IF (fields%name(1:5) == 'insol') THEN ! should may be done in a different way
          ntr_insol = 0
          DO jsp=1,fields%ntr-1
            IF (fields%linsol(jsp))  ntr_insol = ntr_insol + 1
          ENDDO
          IF (ntr_insol == 0) CALL finish('mo_art_init','insoluble tracer in mode '          &
                                &         //TRIM(fields%name)//' is missing.')
          v_conc = pi / 6.0_wp * init_nmb / ntr_insol * fields%info%diameter_ini_nmb**3      &
            &    * fields%info%exp_aero**36 * 1.0E9_wp
          DO jsp=1,fields%ntr-1
            IF (fields%linsol(jsp)) THEN
              fields%info%init_mass_conc(jsp) = v_conc * fields%rho(jsp)
            ELSE
              fields%info%init_mass_conc(jsp) = 0.0_wp
            ENDIF
          ENDDO
        ENDIF
    END SELECT !fields
    current_mode => current_mode%next_mode
  ENDDO !ASSOCIATED(current_mode)

END SUBROUTINE set_meta_init_nmb_mass_conc
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_aerosol_state
