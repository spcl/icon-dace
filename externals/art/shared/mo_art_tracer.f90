!
! mo_art_tracer
! This module provides the tracer definitions. It serves as interface between the reading of XML
! meta data and the definition of the tracers. Definitions of physical tendencies are also
! performed as a part of this module.
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

MODULE mo_art_tracer
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_model_domain,                  ONLY: p_patch
  USE mo_vertical_coord_table,          ONLY: vct_a
  USE mo_exception,                     ONLY: message,finish,message_text
  USE mo_impl_constants,                ONLY: min_rlcell, max_dom, MAX_NTRACER, SUCCESS
  USE mo_var_list,                      ONLY: t_var_list_ptr
  USE mo_var,                           ONLY: t_var
  USE mo_var_metadata_types,            ONLY: t_var_metadata, t_var_metadata_dynamic
  USE mo_fortran_tools,                 ONLY: t_ptr_2d3d
  USE mo_advection_config,              ONLY: t_advection_config
  USE mo_nwp_phy_types,                 ONLY: t_nwp_phy_tend
  USE mo_nonhydro_types,                ONLY: t_nh_prog
  USE mo_art_config,                    ONLY: art_config,nart_tendphy, &
                                          &   ctracer_art
  USE mo_util_string,                   ONLY: split_string
  USE mo_nonhydrostatic_config,         ONLY: kstart_tracer,htop_aero_proc
  USE mo_key_value_store,               ONLY: t_key_value_store
! ART
  USE mo_art_aerosol_state,             ONLY: art_aerosol_state, art_assign_tracers_to_modes, &
                                          &   set_meta_init_nmb_mass_conc
  USE mo_art_modes_linked_list,         ONLY: p_mode_state, t_mode
  USE mo_art_data,                      ONLY: p_art_data, t_art_aero
  USE mo_art_atmo_data,                 ONLY: t_art_atmo
  USE mo_art_impl_constants,            ONLY: IART_AERO_TR, IART_CHEM_TR,   &
                                          &   IART_VARNAMELEN, IART_XMLTAGLEN
  USE mo_art_tracer_def_wrapper,        ONLY: art_tracer_def_wrapper
  USE mo_art_read_xml,                  ONLY: t_xml_file, art_open_xml_file, art_close_xml_file,  &
                                          &   art_get_childnumber_xml, art_read_elements_xml
  USE mo_art_tagging,                   ONLY: get_number_tagged_tracer
  USE mo_art_chem_init_meta,            ONLY: art_chem_init_meta
  USE mo_art_setup_chem_productions,    ONLY: art_setup_chem_productions
  USE mo_art_init_chemtracer,           ONLY: art_get_chemtracer_idx
  USE mo_art_emiss_state,               ONLY: art_aerosol_assign_emission2tracer
  USE mo_art_string_tools,              ONLY: key_value_storage_as_string

#include "add_var_acc_macro.inc"

  IMPLICIT NONE
  
  PRIVATE
  
  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_tracer'
  
  PUBLIC  :: art_tracer
  
  CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_tracer(defcase,jg,nblks_c,this_list,vname_prefix,                &
  &                   ptr_arr,advconf,phy_tend,p_prog,timelev,ldims)
!<
! SUBROUTINE art_tracer
! This subroutine performs the tracer definitions. It serves as interface between reading of XML
! meta data, the definition of the tracers and the definitions of physical tendencies.
! Part of Module: mo_art_tracer
! Author: Daniel Rieger, KIT
! Initial Release: 2016-08-01
! Modifications:
! YYYY-MM-DD: <name>,<institution>
! - <description>
!>
  CHARACTER(LEN=*), INTENT(in)                :: &
    &  defcase                                     !< cases:  prog, conv, turb
  INTEGER, INTENT(in)                         :: &
    &  jg,                                       & !< patch id
    &  nblks_c                                     !< number of blocks (dimension)
  TYPE(t_var_list_ptr), INTENT(inout)             :: &
    &  this_list                                   !< current list  
  CHARACTER(LEN=*), INTENT(in)                :: &  
    &  vname_prefix                                !< prefix for variable names (usually none)
  TYPE(t_ptr_2d3d), INTENT(inout)             :: &
    &  ptr_arr(:)                                  !<
  TYPE(t_advection_config),INTENT(inout)      :: &
    &  advconf                                     !< advection configuration
  TYPE(t_nwp_phy_tend),INTENT(inout),OPTIONAL :: &
    &  phy_tend                                    !<
  TYPE(t_nh_prog), INTENT(inout), OPTIONAL    :: &
    &  p_prog                                      !<
  INTEGER, INTENT(in), OPTIONAL               :: &
    &  timelev,                                  &  !< Only necessary for prognostic tracer list.
    &  ldims(3)                                     !< local dimensions, for checking
! Local Variables
  TYPE(t_var_metadata_dynamic), POINTER   :: &
    &  info_dyn                                !< returns reference to dynamic tracer metadata 
                                               !  of current element
  TYPE(t_xml_file)                        :: &
    &  chem_xmlfile, aero_xmlfile              !< Chemistry/aerosol tracer XML file
                                               !  contains path and file handle
  TYPE(t_key_value_store)                 :: &
    &  storage                                 !< Key-Value storage to be filled with metadata
  INTEGER                                 :: &
    &  jt,                                   & !< counter for number of tagged tracer
    &  jk,jk1,                               & !< loop indice for vertical loop
    &  tracer_idx_icon,                      & !< Index of tracer in ICON this_list container
    &  ntracer_xml,                          & !< Number of tracer in XML file
    &  idx_tracer_xml,                       & !< Index of tracer in XML file
    &  ierror,                               & !< error index of storage
    &  iv,                                   & !< loop index
    &  nmodes,                               & !< number of modes in <mode>-Tag of aerosol-XML
    &  imodes,                               & !<  ^corresponding loop variable
    &  label_i,                              & !< INDEX(labels,",")
    &  modeNumber_i                            !< INDEX(modeNumbers,",")
  LOGICAL                                 :: &
    &  label_given                             !< is label given in XML
  CHARACTER(:), ALLOCATABLE               :: &
    &  tracer_name,                          & !< Name of tracer
    &  tag_name,                             & !< Name of tag to be added to tracer name
    &  c_tmp
  CHARACTER(LEN=IART_VARNAMELEN)          :: &
    &  tracer_name_mapping                     !< Name of tracer for mapping routine (fixed length 
  CHARACTER(LEN=2)                        :: &
    &  cjg                                     !< patch id as character
                                               !  as cce bugfix
  CHARACTER(LEN=100)                      :: &
    &  tracer_name_in,                       &
    &  label_in,                             &
    &  modeNumber_in
  CHARACTER(LEN=3)                        :: &
    &  idx_tracer_str, jt_string               !< idx_tracer_xml / jt as character
  INTEGER, ALLOCATABLE ::                    &
    &  starts(:),                            & !< start-indices of individual modes in modestring
    &  lens(:)                                 !< lengths of individual modes in modestring
  CHARACTER(len=IART_XMLTAGLEN)           :: &
    &  mode,                                 & !< (list of) Names of the modes the aerosol tracer is
    &  labels,                               & !   contained in and the corresponding labels
    &  modeNumbers                             !                and the corresponding modeNumbers (GRIB2)
  CHARACTER(len=IART_VARNAMELEN)          :: &
    &  mode_in

  REAL(wp)                                :: &
    &  htop_tracer_proc                        !< Top height (in m) of the model domain where this
                                               !< individual ART tracer is 
                                               !< transported/diffused/modified
                                               !< (=htop_proc in tracer.xml file)
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo
  TYPE(t_art_aero), POINTER :: &
    &  art_aero
! ----------------------------------
! --- 1.0 Initializations
! ----------------------------------

  ! Set number of convection/turbulent tracers to zero
  ! (this is necessary as art_tracer is called twice, once per timelevel 
  !   and the number of conv/turb would then count double)

  art_atmo => p_art_data(jg)%atmo
  art_aero => p_art_data(jg)%aero

  IF(TRIM(defcase) /= 'conv') art_config(jg)%nturb_tracer = 0
  IF(TRIM(defcase) /= 'turb') art_config(jg)%nconv_tracer = 0

  IF(PRESENT(timelev) .AND. TRIM(defcase) == 'prog') THEN
    IF (timelev == 1) CALL p_art_data(jg)%dict_tracer%init(.FALSE.)
  ENDIF

! ----------------------------------
! --- 2.0 Read tracer metadata and define tracers
! ----------------------------------

  ! ----------------------------------
  ! --- 2.0.1 Open emission XML
  ! ----------------------------------

! ----------------------------------
! --- 2.1 Aerosol tracers
! ----------------------------------

  IF (art_config(jg)%lart_aerosol) THEN
    CALL art_open_xml_file(TRIM(art_config(jg)%cart_aerosol_xml),aero_xmlfile)

    CALL art_get_childnumber_xml(aero_xmlfile, '/tracers/', ntracer_xml)

    DO idx_tracer_xml = 1, ntracer_xml
      ! Create a storage container
      CALL storage%init(.FALSE.)

      ! Read name and metadata of current tracer in XML file
      WRITE(idx_tracer_str,'(I3)') idx_tracer_xml
      CALL art_read_elements_xml(aero_xmlfile,'/tracers/*['//TRIM(ADJUSTL(idx_tracer_str))//']/', &
        &                        storage)
      CALL key_value_storage_as_string(storage,'name',tracer_name,ierror)

      IF(ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_tracer','Tracer name not found.')

      !Get mode information to combine final tracer names
      CALL key_value_storage_as_string(storage,'mode',c_tmp,ierror)
      IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_tracer_def_wrapper',                &
                               &  'Required aerosol metadata mode not present:'//TRIM(tracer_name))
      WRITE(mode,'(A)') TRIM(c_tmp)                         

      ! Get labels information for tracer_mode-combination (optional)
      CALL key_value_storage_as_string(storage,'label',c_tmp,ierror)
      IF (ierror /= SUCCESS) THEN
        CALL message(TRIM(routine)//':art_tracer_def_wrapper',                &
          &          'No label-Tag given for:'//TRIM(tracer_name))
        c_tmp = ""
      END IF
      WRITE(labels,'(A)') TRIM(c_tmp)                         

      ! Get modeNumber information for GRIB2 (optional)
      CALL key_value_storage_as_string(storage,'modeNumber_list',c_tmp,ierror)
      IF (ierror /= SUCCESS) c_tmp = ""
      WRITE(modeNumbers,'(A)') TRIM(c_tmp)                         

      !split modesstring into array of modes
      !maximum number of separated modestrings in string of length N is
      !int((N+1)/2)
      nmodes=FLOOR((LEN_TRIM(mode)+1)/2._wp)
      ALLOCATE(starts(nmodes),lens(nmodes))
      nmodes=1
      CALL split_string(mode,nmodes,starts,lens)
      DO imodes = 1, nmodes
        mode_in = mode(starts(imodes):starts(imodes)+lens(imodes)-1)
        tracer_name_in =  TRIM(tracer_name)//"_"//TRIM(mode_in) ! generic name
        tracer_name_in =  TRIM(ADJUSTL(tracer_name_in))

        !getting label if provided
        label_i = INDEX(labels,',')
        IF(label_i==0) THEN  !no "," found
          IF(LEN(TRIM(labels))>0) THEN ! last entry
            label_in = TRIM(labels)
            labels = ''
          ELSE !no more entries
            label_in = ''
          END IF
        ELSE ! "," found -> get label entry
          label_in = TRIM(labels(1:label_i-1))
        END IF
        
        label_given = .TRUE.
        IF(LEN(TRIM(label_in))==0) THEN ! empty entry -> use tracer_name_in
          label_given = .FALSE.
          label_in = TRIM(tracer_name_in)
        END IF
        labels=labels(label_i+1:LEN(labels)) ! update label-list

        !getting modeNumber for GRIB2 if provided
        modeNumber_i = INDEX(modeNumbers,',')
        IF(modeNumber_i==0) THEN  !no "," found
          IF(LEN(TRIM(modeNumbers))>0) THEN ! last entry
            modeNumber_in = TRIM(modeNumbers)
            modeNumbers = ''
          ELSE !no more entries
            modeNumber_in = ''
          END IF
        ELSE ! "," found -> get modeNumber entry
          modeNumber_in = TRIM(modeNumbers(1:modeNumber_i-1))
        END IF
        
        modeNumbers=modeNumbers(modeNumber_i+1:LEN(modeNumbers)) ! update modeNumber-list

        ! Add tracer and metadata to ICON list
        CALL art_tracer_def_wrapper(IART_AERO_TR,                                       &
          &                         defcase,                                            &
          &                         art_config(jg), advconf,                            &
          &                         this_list,                                          &
          &                         vname_prefix,                                       &
          &                         ptr_arr,                                            &
          &                         storage,                                            &
          &                         aero_xmlfile,                                       &
          &                         '/tracers/*['//TRIM(ADJUSTL(idx_tracer_str))//']/', &
          &                         tracer_idx    = tracer_idx_icon,                    &
          &                         tracer_name_in= label_in,                           &
          &                         tracer_mode   = mode_in,                            &
          &                         tracer_sub    = TRIM(tracer_name),                  &
          &                         modeNumber_in = modeNumber_in,                      &
          &                         timelev       = timelev,                            &
          &                         ldims         = ldims)
  
        IF (TRIM(defcase) == 'prog') THEN
          tracer_name_mapping = TRIM(ADJUSTL(tracer_name_in))
          ctracer_art(tracer_idx_icon) = TRIM(tracer_name_mapping)
          IF(timelev == 1) THEN !MARKER: is this a good strategy?
            CALL p_art_data(jg)%dict_tracer%put(tracer_name_mapping,tracer_idx_icon)
            IF(label_given) CALL p_art_data(jg)%dict_tracer%put(TRIM(ADJUSTL(label_in)),tracer_idx_icon)
          END IF
  
          IF (idx_tracer_xml == 1 .AND. imodes == 1) THEN
            art_aero%itr_start = tracer_idx_icon
            art_aero%itr_end   = tracer_idx_icon + ntracer_xml - 1
          ENDIF
  
          ! modify kstart_tracer according to htop_proc from tracer.xml
          CALL storage%get('htop_proc',htop_tracer_proc,ierror)
          IF (ierror == SUCCESS) THEN
            ! IF htop_proc tag is available in tracer xml
            CALL modify_kstart_tracer(kstart_tracer, jg, tracer_idx_icon, vct_a, htop_tracer_proc)
          ELSE
            ! ELSE use value of htop_aero_proc namelist parameter
            DO jk = 1, art_atmo%nlev
              jk1 = jk + p_patch(jg)%nshift_total
              IF (0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) < htop_aero_proc) THEN
                kstart_tracer(jg,tracer_idx_icon) = jk
                EXIT
              ENDIF
            ENDDO
          ENDIF
          WRITE(message_text,'(a,i4,3a,i4)') 'Domain', jg,         &
            &  '; computations related to ',TRIM(tracer_name_in),    &
            &  ' start in layer ', kstart_tracer(jg,tracer_idx_icon)
          CALL message(routine, message_text)
        ENDIF
      END DO
      ! Checking whether more labels were given than are required
      IF(LEN(TRIM(labels)) > 0) THEN
        WRITE (message_text,*) 'WARNING: More labels provided than modes of species "'//TRIM(tracer_name)// &
          &                    '" in XML-file.'
        CALL message(routine, message_text)
      END IF
      DEALLOCATE(starts,lens)
      ! Set metadata storage free again
      CALL storage%destruct
    ENDDO !idx_tracer_xml

    CALL art_close_xml_file(aero_xmlfile)

    IF (TRIM(defcase) == 'prog') THEN
      IF (timelev == 1) THEN
        write(cjg,"(I2.2)") jg
  
        CALL art_aerosol_state(p_mode_state(jg), TRIM(cjg), TRIM(art_config(jg)%cart_modes_xml), &
          &                    TRIM(art_config(jg)%cart_coag_xml), ldims)
        ! tracer und mode verpointern
        CALL art_assign_tracers_to_modes(p_mode_state(jg),this_list)
        ! write initial number and mass concentration to mode meta data
        CALL set_meta_init_nmb_mass_conc(p_mode_state(jg))

        CALL art_aerosol_assign_emission2tracer(p_mode_state(jg),art_config(jg),this_list,       &
          &                                     p_art_data(jg),                                  &
          &                                     TRIM(art_config(jg)%cart_aero_emiss_xml))
      END IF
    END IF
         
  ENDIF !lart_aerosol

! ----------------------------------
! --- 2.2 Chemical tracers
! ----------------------------------

  IF (art_config(jg)%lart_chem) THEN

    IF (art_config(jg)%lart_chemtracer) THEN
      CALL art_open_xml_file(TRIM(art_config(jg)%cart_chemtracer_xml),chem_xmlfile)
  
      CALL art_get_childnumber_xml(chem_xmlfile, '/tracers/', ntracer_xml)
  
      DO idx_tracer_xml = 1, ntracer_xml
        ! Create a storage container
        CALL storage%init(.FALSE.)
  
        ! Read name and metadata of current tracer in XML file
        WRITE(idx_tracer_str,'(I3)') idx_tracer_xml
        CALL art_read_elements_xml(chem_xmlfile,                                        &
                    &              '/tracers/*['//TRIM(ADJUSTL(idx_tracer_str))//']/',  &
                    &              storage)
  
        CALL key_value_storage_as_string(storage,'name',tracer_name,ierror)
        IF(ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_tracer',    &
                                &         'Tracer name not found.')
  
        IF (get_number_tagged_tracer(storage) /= 0) THEN
          ! There are tagged tracers
          DO jt = 1, get_number_tagged_tracer(storage)
            WRITE(jt_string,'(I3)') jt
            IF (jt < 100) jt_string = '0'//TRIM(ADJUSTL(jt_string))
            IF (jt < 10)  jt_string = '0'//TRIM(ADJUSTL(jt_string))
            CALL key_value_storage_as_string(storage,'tag'//jt_string,tag_name, ierror)
            IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_tracer',    &
                                    &         'tag'//jt_string//' not found.')
            tracer_name_in = TRIM(ADJUSTL(tracer_name))//"_"//TRIM(tag_name)
            ! Add tracer and metadata to ICON list
            CALL art_tracer_def_wrapper(IART_CHEM_TR,                                          &
              &                         defcase,                                               &
              &                         art_config(jg), advconf,                               &
              &                         this_list,                                             &
              &                         vname_prefix,                                          &
              &                         ptr_arr,                                               &
              &                         storage,                                               &
              &                         chem_xmlfile,                                          &
              &                         '/tracers/*['//TRIM(ADJUSTL(idx_tracer_str))//']/',    &
              &                         tracer_idx    = tracer_idx_icon,                       &
              &                         tracer_name_in= tracer_name_in,                        &
              &                         timelev       = timelev,                               &
              &                         ldims         = ldims)
  
            IF (TRIM(defcase) == 'prog') THEN
               ! tracer_name_mapping = TRIM(ADJUSTL(tracer_name_in))       ! in future change to this when all problems caused by this line are solved
               tracer_name_mapping = TRIM(ADJUSTL(tracer_name))
            
               IF(timelev == 1) CALL p_art_data(jg)%dict_tracer%put(tracer_name_mapping,  &
                                                      &             tracer_idx_icon)
  
               ! modify kstart_tracer according to htop_proc from tracer.xml
               CALL storage%get('htop_proc',htop_tracer_proc,ierror)
               IF (ierror == SUCCESS) CALL modify_kstart_tracer(kstart_tracer,     &
                                                     &          jg,                &
                                                     &          tracer_idx_icon,   &
                                                     &          vct_a,             &
                                                     &          htop_tracer_proc)
            ENDIF
          ENDDO
        ELSE
          ! There are no tagged tracer
          tracer_name_in = TRIM(ADJUSTL(tracer_name))
          ! Add tracer and metadata to ICON list
          CALL art_tracer_def_wrapper(IART_CHEM_TR,                                          &
            &                         defcase,                                               &
            &                         art_config(jg), advconf,                               &
            &                         this_list,                                             &
            &                         vname_prefix,                                          &
            &                         ptr_arr,                                               &
            &                         storage,                                               &
            &                         chem_xmlfile,                                          &
            &                         '/tracers/*['//TRIM(ADJUSTL(idx_tracer_str))//']/',    &
            &                         tracer_idx    = tracer_idx_icon,                       &
            &                         tracer_name_in= tracer_name_in,                        &
            &                         timelev       = timelev,                               &
            &                         ldims         = ldims)

          IF (TRIM(defcase) == 'prog') THEN
            tracer_name_mapping = TRIM(ADJUSTL(tracer_name_in))
  
            IF(timelev == 1) CALL p_art_data(jg)%dict_tracer%put(tracer_name_mapping,  &
              &                                                  tracer_idx_icon)
  
            ! modify kstart_tracer according to htop_proc from tracer.xml
            CALL storage%get('htop_proc',htop_tracer_proc,ierror)
            IF (ierror == SUCCESS) CALL modify_kstart_tracer(kstart_tracer,     &
                                                  &          jg,                &
                                                  &          tracer_idx_icon,   &
                                                  &          vct_a,             &
                                                  &          htop_tracer_proc)
          ENDIF
        ENDIF
  
        ! Set metadata storage free again
        CALL storage%destruct
      ENDDO !idx_tracer_xml
  
      CALL art_close_xml_file(chem_xmlfile)
  
    END IF ! lart_chemtracer

    IF (art_config(jg)%lart_mecca) THEN
#ifdef __ART_GPL
      CALL art_open_xml_file(TRIM(art_config(jg)%cart_mecca_xml),chem_xmlfile)
  
      CALL art_get_childnumber_xml(chem_xmlfile, '/tracers/', ntracer_xml)
  
      DO idx_tracer_xml = 1, ntracer_xml
        ! Create a storage container
        CALL storage%init(.FALSE.)
  
        ! Read name and metadata of current tracer in XML file
        WRITE(idx_tracer_str,'(I3)') idx_tracer_xml
        CALL art_read_elements_xml(chem_xmlfile,                                        &
                    &              '/tracers/*['//TRIM(ADJUSTL(idx_tracer_str))//']/',  &
                    &              storage)
  
        CALL key_value_storage_as_string(storage,'name',tracer_name,ierror)
        IF(ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_tracer',    &
                                &         'Tracer name not found.')
  
        IF (get_number_tagged_tracer(storage) /= 0) THEN
          ! There are tagged tracers
          DO jt = 1, get_number_tagged_tracer(storage)
            WRITE(jt_string,'(I3)') jt
            IF (jt < 100) jt_string = '0'//TRIM(ADJUSTL(jt_string))
            IF (jt < 10)  jt_string = '0'//TRIM(ADJUSTL(jt_string))
            CALL key_value_storage_as_string(storage,'tag'//jt_string,tag_name, ierror)
            IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_tracer',    &
                                    &         'tag'//jt_string//' not found.')
  
            tracer_name_in = TRIM(ADJUSTL(tracer_name))//"_"//TRIM(tag_name)
            ! Add tracer and metadata to ICON list
            CALL art_tracer_def_wrapper(IART_CHEM_TR,                                          &
              &                         defcase,                                               &
              &                         art_config(jg), advconf,                               &
              &                         this_list,                                             &
              &                         vname_prefix,                                          &
              &                         ptr_arr,                                               &
              &                         storage,                                               &
              &                         chem_xmlfile,                                          &
              &                         '/tracers/*['//TRIM(ADJUSTL(idx_tracer_str))//']/',    &
              &                         tracer_idx    = tracer_idx_icon,                       &
              &                         tracer_name_in= tracer_name_in,                        &
              &                         timelev       = timelev,                               &
              &                         ldims         = ldims)
  
            IF (TRIM(defcase) == 'prog') THEN
               !tracer_name_mapping = TRIM(ADJUSTL(tracer_name_in))       ! in future change to this when all problems caused by this line are solved
              tracer_name_mapping = TRIM(ADJUSTL(tracer_name))
            
               IF(timelev == 1) CALL p_art_data(jg)%dict_tracer%put(tracer_name_mapping,  &
                                                      &             tracer_idx_icon)
  
               ! modify kstart_tracer according to htop_proc from tracer.xml
               CALL storage%get('htop_proc',htop_tracer_proc,ierror)
               IF (ierror == SUCCESS) CALL modify_kstart_tracer(kstart_tracer,     &
                                                     &          jg,                &
                                                     &          tracer_idx_icon,   &
                                                     &          vct_a,             &
                                                     &          htop_tracer_proc)
            ENDIF
          ENDDO
        ELSE
          ! There are no tagged tracer
          tracer_name_in = TRIM(ADJUSTL(tracer_name))
          ! Add tracer and metadata to ICON list
          CALL art_tracer_def_wrapper(IART_CHEM_TR,                                          &
            &                         defcase,                                               &
            &                         art_config(jg), advconf,                               &
            &                         this_list,                                             &
            &                         vname_prefix,                                          &
            &                         ptr_arr,                                               &
            &                         storage,                                               &
            &                         chem_xmlfile,                                          &
            &                         '/tracers/*['//TRIM(ADJUSTL(idx_tracer_str))//']/',    &
            &                         tracer_idx    = tracer_idx_icon,                       &
            &                         tracer_name_in= tracer_name_in,                        &
            &                         timelev       = timelev,                               &
            &                         ldims         = ldims)

          IF (TRIM(defcase) == 'prog') THEN
            tracer_name_mapping = TRIM(ADJUSTL(tracer_name_in))
  
            IF(timelev == 1) CALL p_art_data(jg)%dict_tracer%put(tracer_name_mapping, &
              &                                                  tracer_idx_icon)
  
            ! modify kstart_tracer according to htop_proc from tracer.xml
            CALL storage%get('htop_proc',htop_tracer_proc,ierror)
            IF (ierror == SUCCESS) CALL modify_kstart_tracer(kstart_tracer,     &
                                                  &          jg,                &
                                                  &          tracer_idx_icon,   &
                                                  &          vct_a,             &
                                                  &          htop_tracer_proc)
          ENDIF
        ENDIF
  
        ! Set metadata storage free again
        CALL storage%destruct
      ENDDO !idx_tracer_xml
  
      CALL art_close_xml_file(chem_xmlfile)
#endif
    END IF ! lart_mecca

    IF (TRIM(defcase) == 'prog') THEN
      CALL art_get_chemtracer_idx(jg)
      CALL art_chem_init_meta(jg,this_list,p_art_data(jg)%chem%param%OH_chem_meta,       &
                   &          art_config(jg)%lart_chemtracer, art_config(jg)%lart_mecca, &
                   &          art_config(jg)%lart_psc)
  
      IF (art_config(jg)%lart_chemtracer) THEN
        CALL art_setup_chem_productions(this_list,p_art_data(jg)%dict_tracer)
      END IF
    END IF
  ENDIF !lart_chem

  !! dumping of dict_tracer, has to be reimplemented in MODULE mo_key_value_store
!  CALL p_art_data(jg)%dict_tracer%dump("p_art_data(jg)%dict_tracer")

! ----------------------------------
! --- 3.0 Collect informations on dimensions
! ----------------------------------

  ! Get the number of tracers that need a physical tendency field (nart_tendphy):
  IF (nart_tendphy == 0) THEN  ! only do this once (This might cause problems 
                               !  when using more than one domain)
    DO iv = 1, this_list%p%nvars
      info_dyn=>this_list%p%vl(iv)%p%info_dyn
      IF (info_dyn%tracer%lis_tracer) THEN
        IF (info_dyn%tracer%lconv_tracer .OR. info_dyn%tracer%lturb_tracer) THEN
          nart_tendphy = nart_tendphy + 1
        ENDIF
      ENDIF
    ENDDO
  ENDIF
  
! ----------------------------------
! --- 4.0 Prepare lists for physical tendency routines
! ----------------------------------
  
  SELECT CASE (TRIM(defcase))
    CASE ('prog')
      WRITE (message_text,*) 'ART: Definition of contiguous tracer fields for convection'  &
        &           //' routine for '//TRIM(this_list%p%vlname)
      CALL message (TRIM(routine)//':art_tracer', message_text)
      CALL define_phy_ptr(jg, art_config(jg)%nconv_tracer, nblks_c, this_list,'conv','prog',    &
        &                   p_prog=p_prog)

      WRITE (message_text,*) 'ART: Definition of contiguous tracer fields for turbulence' &
        &           //' routine for '//TRIM(this_list%p%vlname)
      CALL message (TRIM(routine)//':art_tracer', message_text)
      CALL define_phy_ptr(jg, art_config(jg)%nturb_tracer, nblks_c, this_list,'turb','prog',    &
        &                   p_prog=p_prog)
    CASE ('conv')
      WRITE (message_text,*) 'ART: Definition of contiguous tracer tendencies for convection' &
        &           //' routine for '//TRIM(this_list%p%vlname)
      CALL message (TRIM(routine)//':art_tracer', message_text)
      CALL define_phy_ptr(jg, art_config(jg)%nconv_tracer, nblks_c, this_list,'conv','tend',      &
        &                   phy_tend=phy_tend)
    CASE ('turb')
      WRITE (message_text,*) 'ART: Definition of contiguous tracer tendencies for turbulence'  &
        &           //' routine for '//TRIM(this_list%p%vlname)
      CALL message (TRIM(routine)//':art_tracer', message_text)
      CALL define_phy_ptr(jg, art_config(jg)%nturb_tracer, nblks_c, this_list,'turb','tend',      &
        &                 phy_tend=phy_tend)
    CASE DEFAULT
      CALL finish('mo_art_tracer:art_tracer', &
        &      'defcase unknown')
  END SELECT

! ----------------------------------
! --- 5.0 Cleanup / Success
! ----------------------------------

  CALL message('','ART: All species have been created for defcase '//TRIM(defcase))

  NULLIFY(art_atmo)

END SUBROUTINE art_tracer
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE define_phy_ptr(jg,ntracer_phy,nblks_c,this_list,defcase,ptcase,phy_tend,p_prog)
!<
! SUBROUTINE define_phy_ptr
! For the convection subroutine a continuous field of convective
! tracers has to be passed. Hence, this subroutine creates a list
! of prognostic convective tracers (p_prog%conv_tracer) and tendencies
! of convective tracers (phy_tend%conv_tracer_tend). A relation
! between those convective lists and the prognostic/tendency list
! is stored within p_prog%conv_tracer(:,:)%idx_tracer and
! phy_tend%conv_tracer_tend(:,:)%idx_tracer
! Part of Module: mo_art_tracer
! Author: Daniel Rieger, KIT
! Initial Release: 2012-11-21
! Modifications:
! YYYY-MM-DD: <name>,<institution>
! - <description>
!>
  INTEGER,          INTENT(in)      :: &
    &   jg                                !< patch id
  INTEGER,          INTENT(in)      :: &
    &   ntracer_phy,                   &  !< number of tracers influenced by physical process
    &   nblks_c                           !<
  TYPE(t_var_list_ptr), INTENT(inout)   :: &
    &   this_list                         !< current list
  CHARACTER(LEN=*), INTENT(in)      :: &
    &   defcase                           !< cases: conv, turb
  CHARACTER(LEN=*), INTENT(in)      :: & 
    &   ptcase                            !< cases: prog, tend
  TYPE(t_nwp_phy_tend),INTENT(inout),OPTIONAL :: & 
    &   phy_tend                          !< 
  TYPE(t_nh_prog),     INTENT(inout),OPTIONAL :: &
    &   p_prog                            !< 
! Local variables
  TYPE(t_var_metadata), POINTER     :: &
    &   info                              !< pointer to tracer metadata of current element
  TYPE(t_var_metadata_dynamic), POINTER     :: &
    &   info_dyn                          !< pointer to dynamic tracer metadata of current element
  INTEGER                           :: &
    &    idx_trac, jb,                 &  !< 
    &    i_startblk, i_endblk,         &  !< Index of start and end block 
    &    i_rlstart, i_rlend,           &  !< Start and end values of refined grid
    &    i_nchdom, iv                     !< 
  INTEGER, POINTER                  :: &
    &    jsp                              !< returns index of element


  SELECT CASE(TRIM(defcase))
    CASE ('conv')
      SELECT CASE(TRIM(ptcase))
        CASE ('prog')
          IF (.NOT. PRESENT(p_prog)) THEN
            CALL finish('mo_art_tracer:define_phy_ptr', &
              &      'p_prog not present in ptcase = prog')
          ENDIF
          ALLOCATE (p_prog%conv_tracer(nblks_c,ntracer_phy))
          !$ACC ENTER DATA CREATE(p_prog%conv_tracer)
        CASE ('tend')
          IF (.NOT. PRESENT(phy_tend)) THEN
            CALL finish('mo_art_tracer:define_phy_ptr', &
              &      'phy_tend not present in ptcase = tend')
          ENDIF
          ALLOCATE (phy_tend%conv_tracer_tend(nblks_c,ntracer_phy))
          !$ACC ENTER DATA CREATE(phy_tend%conv_tracer_tend)
      END SELECT
    CASE ('turb')
      SELECT CASE(TRIM(ptcase))
        CASE ('prog')
          IF (.NOT. PRESENT(p_prog)) THEN
            CALL finish('mo_art_tracer:define_phy_ptr', &
              &      'p_prog not present in ptcase = prog')
          ENDIF
          ALLOCATE (p_prog%turb_tracer(nblks_c,ntracer_phy))
          !$ACC ENTER DATA CREATE(p_prog%turb_tracer)
        CASE ('tend')
          IF (.NOT. PRESENT(phy_tend)) THEN
            CALL finish('mo_art_tracer:define_phy_ptr', &
              &      'phy_tend not present in ptcase = tend')
          ENDIF
          ALLOCATE (phy_tend%turb_tracer_tend(nblks_c,ntracer_phy))
          !$ACC ENTER DATA CREATE(phy_tend%turb_tracer_tend)
      END SELECT
    CASE DEFAULT
      CALL finish('mo_art_tracer:define_phy_ptr', &
        &      'Unknown defcase')
  END SELECT
 
  !-----------------------------------
  !>Start the Subroutine
  !-----------------------------------

  i_nchdom   = MAX(1,p_patch(jg)%n_childdom)
  i_rlstart  = 1  !is always one
  i_rlend    = min_rlcell
  i_startblk = p_patch(jg)%cells%start_blk(i_rlstart,1)
  i_endblk   = p_patch(jg)%cells%end_blk(i_rlend,i_nchdom)
  idx_trac   = 1  ! Start index with index 1.
  
  !start DO-loop over elements in list:
  DO iv = 1, this_list%p%nvars
    !get meta data of current element:
    info_dyn=> this_list%p%vl(iv)%p%info_dyn
    info=> this_list%p%vl(iv)%p%info
    jsp => info%ncontained
    
    ! There are four cases
    ! ---- Case 1
    !< If it is a tracer and we want to define the prognostic list
    IF (info_dyn%tracer%lis_tracer) THEN
      IF (info_dyn%tracer%lconv_tracer) THEN !< for convection routine
        IF (TRIM(ptcase) == 'prog' )   THEN
          IF (TRIM(defcase) == 'conv' ) THEN
            DO jb=i_startblk, i_endblk 
              NULLIFY(p_prog%conv_tracer(jb,idx_trac)%ptr)
              p_prog%conv_tracer(jb,idx_trac)%ptr      => this_list%p%vl(iv)%p%r_ptr(:,:,jb,jsp,1)
              __acc_attach(p_prog%conv_tracer(jb,idx_trac)%ptr)
              p_prog%conv_tracer(jb,idx_trac)%idx_tracer =  jsp
            ENDDO
            WRITE (message_text,'(3a,i3,a,i3,a1)') 'ART: Definition for ', info%name,  &
              &                                    'idx_trac = ', idx_trac, ', jsp = ', jsp, '.'
            CALL message('',TRIM(message_text))
            idx_trac=idx_trac+1
          ENDIF
        ENDIF
      ENDIF
    ENDIF
    
    ! ---- Case 2
   !< If it is a tracer and we want to define the prognostic list
    IF (info_dyn%tracer%lis_tracer) THEN
      IF (info_dyn%tracer%lturb_tracer)  THEN !< for turbulence routine
        IF (TRIM(ptcase) == 'prog')   THEN
          IF (TRIM(defcase) == 'turb') THEN
            DO jb=i_startblk, i_endblk 
              NULLIFY(p_prog%turb_tracer(jb,idx_trac)%ptr)
              p_prog%turb_tracer(jb,idx_trac)%ptr       =>this_list%p%vl(iv)%p%r_ptr(:,:,jb,jsp,1)
              __acc_attach(p_prog%turb_tracer(jb,idx_trac)%ptr)
              p_prog%turb_tracer(jb,idx_trac)%idx_tracer= jsp
            ENDDO
            WRITE (message_text,'(3a,i3,a,i3,a1)') 'ART: Definition for ', info%name,           &
              &                                    'idx_trac = ', idx_trac, ', jsp = ', jsp, '.'
            CALL message('',TRIM(message_text))
            idx_trac=idx_trac+1
          ENDIF
        ENDIF
      ENDIF
    ENDIF
    
    ! ---- Case 3
    !< If it is a tendency which is defined as convective tendency
    IF (info_dyn%tracer%lconv_tracer) THEN
      IF (TRIM(ptcase) == 'tend') THEN      !< and if we are at the right cases in the program
        IF (TRIM(defcase) == 'conv') THEN
          DO jb=i_startblk, i_endblk 
            NULLIFY(phy_tend%conv_tracer_tend(jb,idx_trac)%ptr)
            phy_tend%conv_tracer_tend(jb,idx_trac)%ptr  =>this_list%p%vl(iv)%p%r_ptr(:,:,jb,jsp,1)
            __acc_attach(phy_tend%conv_tracer_tend(jb,idx_trac)%ptr)
            phy_tend%conv_tracer_tend(jb,idx_trac)%idx_tracer = jsp
          ENDDO
          WRITE (message_text,'(3a,i3,a,i3,a1)') 'ART: Definition for ', info%name, 'idx_trac = ',&
            &                                    idx_trac, ', jsp = ', jsp, '.'
          CALL message('',TRIM(message_text))
          idx_trac=idx_trac+1
        ENDIF
      ENDIF
    ENDIF
    
    ! ---- Case 4
    !< If it is a tendency which is defined as turbulence tendency
    IF (info_dyn%tracer%lturb_tracer)  THEN
      IF (TRIM(ptcase) == 'tend')   THEN  !< and if we are at the right cases in the program
        IF (TRIM(defcase) == 'turb') THEN
          DO jb=i_startblk, i_endblk 
            NULLIFY(phy_tend%turb_tracer_tend(jb,idx_trac)%ptr)
            phy_tend%turb_tracer_tend(jb,idx_trac)%ptr  =>this_list%p%vl(iv)%p%r_ptr(:,:,jb,jsp,1)
            __acc_attach(phy_tend%turb_tracer_tend(jb,idx_trac)%ptr)
            phy_tend%turb_tracer_tend(jb,idx_trac)%idx_tracer = jsp
          ENDDO
          WRITE(message_text,'(3a,i3,a,i3,a1)') 'ART: Definition for ', info%name, 'idx_trac = ', &
            &                                    idx_trac, ', jsp = ', jsp, '.'
          CALL message('',TRIM(message_text))
          idx_trac=idx_trac+1
        ENDIF
      ENDIF
    ENDIF

  ENDDO !loop iv

END SUBROUTINE define_phy_ptr
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE modify_kstart_tracer(kstart_tracer, jg, tracer_idx_icon, vct_a, htop_tracer_proc)
  !<
  ! SUBROUTINE modify_kstart_tracer
  ! Computes kstart_tracer, the top level for the individual ART-tracer
  ! transport/diffusion/modification
  !
  ! Initially, kstart_tracer(jg,MAX_NTRACER) holds the same kstart level for all tracers;
  ! this value is derrived from htop_tracer_proc, which may be set in the nonhydrostatic_nml.
  ! If htop_proc is specified as xml-tag in the tracers.xml file for a specific tracer, this
  ! subroutine is called and overwrites the kstart_tracer value with the specific level derived
  ! from htop_proc
  !
  ! Based on: -
  ! Part of Module: mo_art_tracer
  ! Author: Andrea Steiner, DWD
  ! Initial Release: 2018-08-17
  ! Modifications:
  ! 2018-08-17: Andrea Steiner, DWD
  ! -
  !>
  INTEGER, INTENT(inout)                      :: &
    &  kstart_tracer(max_dom,MAX_NTRACER)          !< start level for art tracers
  INTEGER, INTENT(in)                         :: &
    &  jg,                                       & !< patch id
    &  tracer_idx_icon                             !< Index of tracer in ICON this_list container
  REAL(wp)                                    :: &
    &  vct_a(:),                                 & !<param. A of the vertical coordinate
    &  htop_tracer_proc                            !< Top height (in m) of the model domain where 
                                                   !< this individual ART tracer is 
                                                   !< transported/diffused/modified
                                                   !< (=htop_proc in tracer.xml file)
  ! Local Variables
  INTEGER                                     :: &
    &  nlev,                                     & !< number of vertical levels
    &  nshift_total,                             & !< total shift of model top with respect 
    &  jk,jk1                                      !< to global domain loop indices


  ! number of vertical levels
  nlev         = p_patch(jg)%nlev
  nshift_total = p_patch(jg)%nshift_total

  DO jk = 1, nlev
     jk1 = jk + nshift_total
     IF (0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) < htop_tracer_proc) THEN
        kstart_tracer(jg,tracer_idx_icon) = jk
        EXIT
     ENDIF
  ENDDO

END SUBROUTINE modify_kstart_tracer
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_tracer


