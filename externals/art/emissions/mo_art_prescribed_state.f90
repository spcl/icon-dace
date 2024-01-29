!
! mo_art_prescribed_state
! This module provides data structures for general prescribed data
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

MODULE mo_art_prescribed_state
  ! ICON
  USE mo_kind,                     ONLY: wp
  USE mo_mpi,                      ONLY: p_max, p_min, p_comm_work
  USE mo_exception,                ONLY: finish, message, message_text
  USE mo_var_list,                 ONLY: t_var_list_ptr
  USE mo_var,                      ONLY: t_var
  USE mo_impl_constants,           ONLY: SUCCESS
  USE mo_art_config,               ONLY: IART_PATH_LEN, art_config
  USE mo_tracer_metadata_types,    ONLY: t_chem_meta, t_aero_meta
  USE mo_key_value_store,          ONLY: t_key_value_store
  USE mo_model_domain,             ONLY: t_patch
  USE mo_time_config,              ONLY: time_config
  USE mo_util_vgrid_types,         ONLY: vgrid_buffer
  ! EXTERNALS
  USE mtime,                       ONLY: max_datetime_str_len, newDatetime, &
                                     &   datetime
  ! ART
  USE mo_art_data,                 ONLY: p_art_data
  USE mo_art_atmo_data,            ONLY: t_art_atmo
  USE mo_art_wrapper_routines,     ONLY: art_get_indices_c
  USE mo_art_prescribed_types,     ONLY: t_art_prescr_list,         &
                                     &   t_art_prescr_list_element, &
                                     &   t_art_prescr_var_dep
  USE mo_art_read_emissions,       ONLY: art_init_emission_struct,  &
                                     &   art_read_emissions
  USE mo_art_read_xml,             ONLY: t_xml_file,          &
                                     &   art_get_prescr_meta, &
                                     &   art_open_xml_file,   &
                                     &   art_close_xml_file
  USE mo_art_impl_constants,       ONLY: IART_VARNAMELEN
  USE mo_art_io_constants,         ONLY: IEXT_CHEM_SPEC
  USE mo_art_create_filenames,     ONLY: art_get_res_string
  USE mo_art_chem_init_utils,      ONLY: deallocate_chem_init_chem, &
                                     &   deallocate_chem_init_atm
  USE mo_art_string_tools,         ONLY: key_value_storage_as_string
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_create_prescr_list
  PUBLIC :: art_prescribe_tracers
  PUBLIC :: art_delete_prescr_list
CONTAINS

SUBROUTINE art_create_prescr_list(this_list,p_prog_list,jg)
!<
! SUBROUTINE art_create_prescr_list
! This subroutine creates the linked list with meta data 
! of the dataset for prescribing tracers by reading the meta data from the
! key_value_store and the XML file for external datasets
! Part of Module: mo_art_prescribed_state
! Author: Michael Weimer, KIT
! Initial Release: 2017-09-06
! Modifications:
!>
  IMPLICIT NONE
  TYPE(t_art_prescr_list), INTENT(INOUT), TARGET :: &
    &  this_list                 !< linked list to be created
  TYPE(t_var_list_ptr), INTENT(IN) :: &
    &  p_prog_list               !< list of prognostic tracers
  INTEGER, INTENT(IN)       :: &
    &  jg                        !< patch on which computation is performed
  ! local variables
  TYPE(t_art_prescr_list_element), POINTER :: &
    &  prescr_list_element,  &   !< element of this_list
    &  prescr_list_element2      !< element of this_list
  CHARACTER(:), ALLOCATABLE :: &
    &  dataset_name,                  & !< name of the dataset
    &  var_name_in_dataset              !< variable name in the dataset to be
                                        !  connected with the tracer to be prescribed
  CHARACTER(LEN = IART_VARNAMELEN) :: &
    &  dataset_idx_str                  !< character of dataset_idx
  INTEGER :: &
    &  ierror_dataset,     &     !< error index of the key_value_store, loop indices
    &  ierror, jsp, i,     &
    &  dataset_idx, j, iv
  LOGICAL :: &
    &  dataset_found             !< flag if dataset is found
  TYPE(t_xml_file) :: &
    &  tixi_file                 !< XML file containing the "global" meta data of the dataset
  TYPE(t_key_value_store), POINTER ::  &
    &  meta_storage              !< pointer to the tracer storage
  TYPE(datetime), POINTER :: &
    &  current_mtime             !< current simulation time
  TYPE(t_art_prescr_var_dep) :: &
    &  var_dep_save              !< saved values of variable dependent meta information
  INTEGER ::             &
    &  jb, jc, nlen,     &      !< loop indices
    &  minloc_zifc(3),   &      !< location in array of minimum of z_ifc
    &  bot_idx,          &      !< upper and bottom indices of the prescribing
    &  upp_idx,          &
    &  bot2, upp2
  REAL(wp), ALLOCATABLE :: &
    &  zifc_sea(:)              !< actual minimum value of z_ifc
  REAL(wp)                      ::   &
    &  mol_weight,                   &  !< molar weight of species (kg / mol)
    &  bottom_height,                &  !< bottom height and upper height of the prescription
    &  upper_height
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo                !< pointer to ART atmo fields

  TYPE(t_chem_meta),POINTER :: chem_meta_ptr
  TYPE(t_aero_meta),POINTER :: aero_meta_ptr

  art_atmo => p_art_data(jg)%atmo

  ALLOCATE(zifc_sea(art_atmo%nlevp1))

  this_list%first_prescr_element => NULL()
  this_list%num_elements = 0

  IF (TRIM(art_config(jg)%cart_ext_data_xml) /= '') THEN
    current_mtime => time_config%tc_current_date

    CALL art_open_xml_file(art_config(jg)%cart_ext_data_xml,tixi_file)

    ! go through all tracer elements and look into their storages if tracers
    ! should be prescribed
    DO iv = 1, p_prog_list%p%nvars
      dataset_found = .FALSE.

      IF (p_prog_list%p%vl(iv)%p%info_dyn%tracer%lis_tracer) THEN
        jsp = p_prog_list%p%vl(iv)%p%info%ncontained

        ! set pointer to the tracer storage or skip the tracers with wrong type
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

        ! read the dataset name and variable name
        CALL key_value_storage_as_string(meta_storage,'dataset01',dataset_name,ierror_dataset)
        dataset_idx = 1
        WRITE(dataset_idx_str,'(I2.2)') dataset_idx

        DO WHILE (ierror_dataset == SUCCESS)

          CALL meta_storage%get('mol_weight',mol_weight,ierror)
          IF (ierror /= SUCCESS) THEN
            CALL finish('mo_art_prescribed_state:art_create_prescr_list',          &
                    &   'For prescribing '//TRIM(p_prog_list%p%vl(iv)%p%info%name)  &
                    & //', mol_weight has to be available.')
          END IF

          CALL key_value_storage_as_string(meta_storage,                     &
                     &          'cvar_name@dataset'//TRIM(dataset_idx_str),  &
                     &          var_name_in_dataset,                         &
                     &          ierror)

          IF (ierror == SUCCESS) THEN

            minloc_zifc = (/ 1,art_atmo%nlev+1,1 /)
            DO jb = art_atmo%i_startblk,art_atmo%i_endblk
              IF (jb < art_atmo%nblks) THEN
                nlen = art_atmo%nproma
              ELSE
                nlen = art_atmo%npromz
              END IF

              DO jc = 1,nlen
                IF (art_atmo%z_ifc(jc,art_atmo%nlevp1,jb)  &
                   & < art_atmo%z_ifc(minloc_zifc(1),minloc_zifc(2),minloc_zifc(3))) THEN

                  minloc_zifc = (/ jc,art_atmo%nlevp1,jb /)
                END IF
              END DO
            END DO

            zifc_sea(:) = art_atmo%z_ifc(minloc_zifc(1),:,minloc_zifc(3))

            CALL meta_storage%get('rbottom_height@dataset'//TRIM(dataset_idx_str),  &
                      &           bottom_height,                                    &
                      &           ierror)

            IF (ierror /= SUCCESS) THEN
              CALL message('mo_art_prescribed_state:art_create_prescr_list',              &
                     &    'WARNING: attribute rbottom_height missing for dataset with '   &
                     &  //'variable '//TRIM(var_name_in_dataset)//' for tracer '          &
                     &  //TRIM(p_prog_list%p%vl(iv)%p%info%name)//'. Set to surface (0 m).')
              bottom_height = 0._wp
            END IF

            CALL meta_storage%get('rupper_height@dataset'//TRIM(dataset_idx_str),  &
                       &          upper_height,                                    &
                       &          ierror)

            IF (ierror /= SUCCESS) THEN
              CALL message('mo_art_prescribed_state:art_create_prescr_list',              &
                     &    'WARNING: attribute rupper_height missing for dataset with '    &
                     &  //'variable '//TRIM(var_name_in_dataset)//' for tracer '          &
                     &  //TRIM(p_prog_list%p%vl(iv)%p%info%name)//'. Set to top of atmosphere.')
              upper_height = MAXVAL(zifc_sea)
            END IF

            CALL art_set_bot_upp_indices(zifc_sea,upper_height,bottom_height,  &
                       &                 bot_idx,upp_idx)

            ! first look into this_list if it already exists
            IF (.NOT. ASSOCIATED(this_list%first_prescr_element)) THEN
              ALLOCATE(this_list%first_prescr_element)
              prescr_list_element => this_list%first_prescr_element
            ELSE
              prescr_list_element => this_list%first_prescr_element

              DO i = 1,this_list%num_elements
                IF (TRIM(prescr_list_element%prescr%id) == TRIM(dataset_name)) THEN
                  dataset_found = .TRUE.
                  EXIT
                END IF
      
                IF (i < this_list%num_elements) THEN
                  prescr_list_element => prescr_list_element%next_prescr_element
                END IF
              END DO

              IF (.NOT. dataset_found) THEN
                ALLOCATE(prescr_list_element%next_prescr_element)
                prescr_list_element => prescr_list_element%next_prescr_element
                prescr_list_element%next_prescr_element => NULL()
              END IF
            END IF

            prescr_list_element%prescr%prescribed_state%chem_init_chem%n_spec_readin  &
              &   =>  prescr_list_element%prescr%num_vars

            IF (dataset_found) THEN
              ! if the dataset is already in the list, just put the new variable
              ! name and the tracer index
              CALL art_add_var_to_dataset(prescr_list_element,    &
                       &                  var_name_in_dataset,    &
                       &                  jsp, bot_idx, upp_idx,  &
                       &                  mol_weight)
            ELSE
              ! if not then look into xml if dataset is found
              CALL art_new_prescr_element(prescr_list_element,jg,       &
                        &                 tixi_file,dataset_name,       &
                        &                 var_name_in_dataset,          &
                        &                 jsp, bot_idx, upp_idx,        &
                        &                 mol_weight)

              this_list%num_elements = this_list%num_elements + 1
            END IF
          ELSE
            CALL finish('mo_art_prescribed_state:art_create_prescr_list',  &
                   &    'could not find cvar_name attribute in dataset '   &
                   &  //'element for '//TRIM(p_prog_list%p%vl(iv)%p%info%name))
          END IF
        
          dataset_idx = dataset_idx + 1
          WRITE(dataset_idx_str,'(I2.2)') dataset_idx
          CALL key_value_storage_as_string(meta_storage,'dataset'//TRIM(dataset_idx_str),dataset_name,ierror_dataset)
        END DO
      END IF

    END DO

    CALL art_close_xml_file(tixi_file)

    NULLIFY(current_mtime)
  END IF  ! cart_ext_data_xml /= ''

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Error checking of for height ranges
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  prescr_list_element => this_list%first_prescr_element
  DO WHILE (ASSOCIATED(prescr_list_element))

    DO i = 1,prescr_list_element%prescr%num_vars
      var_dep_save = prescr_list_element%prescr%var_dep(i)

      IF (i+1 <= prescr_list_element%prescr%num_vars) THEN
        DO j = i+1,prescr_list_element%prescr%num_vars
          IF (var_dep_save%tracer_idx == prescr_list_element%prescr%var_dep(j)%tracer_idx) THEN
            bot2 = prescr_list_element2%prescr%var_dep(j)%bot_idx
            upp2 = prescr_list_element2%prescr%var_dep(j)%upp_idx

            IF (((var_dep_save%bot_idx >= upp2) .AND. (var_dep_save%bot_idx <= bot2))  &
             &  .OR. ((var_dep_save%upp_idx >= upp2) .AND. (var_dep_save%upp_idx <= bot2))) THEN
               WRITE(message_text,'(A,I3,A)') &
                      &    'The given heights for dataset with variable name ' &
                      &   //TRIM(prescr_list_element%prescr%vname(j))          &
                      &   //' overlap with another for tracer with id ',       &
                      &   var_dep_save%tracer_idx, '.'
               CALL finish('mo_art_prescribed_state:art_create_prescr_list', message_text)
            END IF ! error check
          END IF ! tracer indices equal?
        END DO ! number of variables (j)
      END IF ! If do loop is possible

      prescr_list_element2 => prescr_list_element%next_prescr_element

      DO WHILE (ASSOCIATED(prescr_list_element2))
        DO j = 1,prescr_list_element2%prescr%num_vars
          IF (var_dep_save%tracer_idx == prescr_list_element2%prescr%var_dep(j)%tracer_idx) THEN
            bot2 = prescr_list_element2%prescr%var_dep(j)%bot_idx
            upp2 = prescr_list_element2%prescr%var_dep(j)%upp_idx

            IF (((var_dep_save%bot_idx >= upp2) .AND. (var_dep_save%bot_idx <= bot2))  &
             &  .OR. ((var_dep_save%upp_idx >= upp2) .AND. (var_dep_save%upp_idx <= bot2))) THEN
               WRITE(message_text,'(A,I3,A)') &
                      &    'The given heights for dataset with variable name ' &
                      &   //TRIM(prescr_list_element%prescr%vname(j))          &
                      &   //' overlap with another for tracer with id ',       &
                      &   var_dep_save%tracer_idx, '.'
               CALL finish('mo_art_prescribed_state:art_create_prescr_list', message_text)
            END IF ! error check
          END IF ! tracer indices equal?
        END DO ! number of variables (j)
        prescr_list_element2 => prescr_list_element2%next_prescr_element
      END DO  ! Associated (element2)

    END DO ! number of variables (i)
    prescr_list_element => prescr_list_element%next_prescr_element
  END DO ! associated (element)

  DEALLOCATE(zifc_sea)

END SUBROUTINE art_create_prescr_list
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_add_var_to_dataset(prescr_list_element,var_name_in_dataset,  &
                        &         tracer_idx, bot_idx, upp_idx, mol_weight)
!<
! SUBROUTINE art_add_var_to_dataset
! This subroutine adds a variable to an existing element of the list for
! prescribing tracers
! Part of Module: mo_art_prescribed_state
! Author: Michael Weimer, KIT
! Initial Release: 2017-09-06
! Modifications:
!>
  IMPLICIT NONE
  TYPE(t_art_prescr_list_element), POINTER :: &
    &  prescr_list_element        !< element in linked list containing the meta data of the dataset
  CHARACTER(LEN=*), INTENT(IN) :: &
    &  var_name_in_dataset        !< variable name in the dataset
  INTEGER, INTENT(IN) :: &
    &  tracer_idx,       &        !< index of the tracer
    &  bot_idx,          &        !< bottom index
    &  upp_idx                    !< upper index
  REAL(wp), INTENT(IN) :: &
    &  mol_weight                 !< molar weight in kg / mol
  ! local variable
  CHARACTER(LEN=IART_VARNAMELEN), ALLOCATABLE :: &
    &  var_names_saved(:)         !< temporal saves of the current arrays
  TYPE(t_art_prescr_var_dep), ALLOCATABLE :: &       !  (they are increased by one and the parameter
    &  var_dep_saved(:)    !   are appended to the increased array)

  ! write current arrays into the "saved" arrays and deallocate the old arrays
  ALLOCATE(var_names_saved(prescr_list_element%prescr%num_vars))
  var_names_saved = prescr_list_element%prescr%vname(:)
  DEALLOCATE(prescr_list_element%prescr%vname)

  ALLOCATE(var_dep_saved(prescr_list_element%prescr%num_vars))
  var_dep_saved = prescr_list_element%prescr%var_dep(:)
  DEALLOCATE(prescr_list_element%prescr%var_dep)


  ! increase number of variables by 1
  prescr_list_element%prescr%num_vars = prescr_list_element%prescr%num_vars + 1


  ! reallocate the arrays and put the old values in 1:num_vars-1 and new value
  ! into the last element of the array
  ALLOCATE(prescr_list_element%prescr%var_dep(  &
                &         prescr_list_element%prescr%num_vars))
  prescr_list_element%prescr%var_dep(1:prescr_list_element%prescr%num_vars - 1)  &
                            &       = var_dep_saved(:)
  prescr_list_element%prescr%var_dep(prescr_list_element%prescr%num_vars)%tracer_idx &
                            &       = tracer_idx
  prescr_list_element%prescr%var_dep(prescr_list_element%prescr%num_vars)%bot_idx &
                            &       = bot_idx
  prescr_list_element%prescr%var_dep(prescr_list_element%prescr%num_vars)%upp_idx &
                            &       = upp_idx
  prescr_list_element%prescr%var_dep(prescr_list_element%prescr%num_vars)%mol_weight &
                            &       = mol_weight

  ALLOCATE(prescr_list_element%prescr%vname(prescr_list_element%prescr%num_vars))
  prescr_list_element%prescr%vname(1:prescr_list_element%prescr%num_vars - 1)  &
                     &              = var_names_saved(:)

  prescr_list_element%prescr%vname(prescr_list_element%prescr%num_vars)   &
                     &              = TRIM(var_name_in_dataset)


  DEALLOCATE(var_dep_saved)
  DEALLOCATE(var_names_saved)
END SUBROUTINE art_add_var_to_dataset
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_new_prescr_element(prescr_list_element,jg,tixi_file,   &
                         &        dataset_name,var_name_in_dataset,   &
                         &        tracer_idx, bot_idx, upp_idx,       &
                         &        mol_weight)
!<
! SUBROUTINE art_new_prescr_element
! This subroutine creates a new element in the linked list with meta data for
! prescribing tracers. An element contains global information such as the name
! of the dataset as well as specific information about the variables in this
! dataset which will be used for prescribing
! Part of Module: mo_art_prescribed_state
! Author: Michael Weimer, KIT
! Initial Release: 2017-09-06
! Modifications:
!>
  IMPLICIT NONE
  TYPE(t_art_prescr_list_element), POINTER :: &
    &  prescr_list_element     !< element of the linked list to be created
  INTEGER, INTENT(IN) :: &
    &  jg                      !< p_patch%id
  TYPE(t_xml_file), INTENT(IN) :: &
    &  tixi_file               !< XML file with global meta info of the dataset
  CHARACTER(LEN=*), INTENT(IN) :: &
    &  dataset_name,              & !< name of the dataset
    &  var_name_in_dataset          !< name of the variables in the dataset
  INTEGER, INTENT(IN) :: &
    &  tracer_idx,       &       !< index of the tracer
    &  bot_idx,          &       !< indices of the prescription
    &  upp_idx
  REAL(wp), INTENT(IN) :: &
    &  mol_weight                !< molar weight in kg / mol
  ! local variables
  CHARACTER(LEN = max_datetime_str_len) :: &
    &  first_date, last_date   !< first and last date of the dataset
  CHARACTER(LEN = IART_PATH_LEN) :: &
    &  vert_coord              !< name of the coordinate file
  CHARACTER(LEN = IART_VARNAMELEN) :: &
    &   vert_unit,                    & !< unit of vertical coordinate
    &   model_name                      !< name of the model
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo

  art_atmo => p_art_data(jg)%atmo

  ! read the meta data from the xml file
  CALL art_get_prescr_meta(tixi_file,dataset_name,var_name_in_dataset, &
                   &       first_date,last_date,vert_coord,vert_unit,model_name)

  ! initialise the list element
  prescr_list_element%prescr%num_vars = 1
  ALLOCATE(prescr_list_element%prescr%vname(1),      &
           prescr_list_element%prescr%var_dep(1))

  ! write the data from XML to the new element
  prescr_list_element%prescr%var_dep(1)%tracer_idx = tracer_idx
  prescr_list_element%prescr%var_dep(1)%bot_idx = bot_idx
  prescr_list_element%prescr%var_dep(1)%upp_idx = upp_idx
  prescr_list_element%prescr%var_dep(1)%mol_weight = mol_weight
  prescr_list_element%prescr%vname(1) = TRIM(var_name_in_dataset)
  prescr_list_element%prescr%first_date => newDatetime(first_date)
  prescr_list_element%prescr%last_date => newDatetime(last_date)
  prescr_list_element%prescr%unit_vertical = TRIM(vert_unit)
  prescr_list_element%prescr%iType_data = IEXT_CHEM_SPEC
  prescr_list_element%prescr%vert_coord_filename = vert_coord
  prescr_list_element%prescr%id = TRIM(dataset_name)

  prescr_list_element%prescr%prescribed_state%chem_init_in%model_name   &
                        &       = TRIM(model_name)

  prescr_list_element%prescr%path = TRIM(art_config(jg)%cart_input_folder)  &
    &    //'/datasets/'//TRIM(dataset_name)//'/'                            &
    &    //TRIM(art_get_res_string(jg))                                     &
    &    //'_'//TRIM(art_config(jg)%cart_io_suffix)

  NULLIFY(art_atmo)
  
END SUBROUTINE art_new_prescr_element
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_set_bot_upp_indices(zifc_sea,upper_height,bottom_height,bot_idx,upp_idx)
!<
! SUBROUTINE art_set_bot_upp_indices
! This subroutine computes the lower and upper levels within which the tracer
! should be prescribed. It uses z_ifc over a sea point and searches until the
! XML given upper or lower height is exceeded. Then, the height is between this
! and next upper or lower half level. The index is chosen so that the full
! levels between these half level will be chosen as prescribed height range.
! Part of Module: mo_art_prescribed_state
! Author: Michael Weimer, KIT
! Initial Release: 2017-12-06
! Modifications:
!>
  IMPLICIT NONE
  REAL(wp), INTENT(IN) :: &
    &  zifc_sea(:)          !< model interface level height above a grid point over sea (m)
  REAL(wp), INTENT(IN) :: &
    &  upper_height,      & !< given upper height bound below which tracer is prescribed (m)
    &  bottom_height        !< given lower height bound above which tracer is prescribed (m)
  INTEGER, INTENT(INOUT) :: &
    &  bot_idx,           & !< (full) level index corresponding to bottom_height
    &  upp_idx              !< (full) level index corresponding to upper_height 
                            !  (lower than bot_idx!!)
  ! local variables
  INTEGER ::        &
    &  nlev,        &       !< number of full levels
    &  bot_idx_max, &       !< global maximum bottom height index
    &  upp_idx_min          !< global minimum upper height index


  nlev = SIZE(zifc_sea) - 1
  bot_idx = nlev
  upp_idx = 2

  DO WHILE ((bot_idx > 1)   &
       &     .AND. (zifc_sea(bot_idx) < bottom_height))
    bot_idx = bot_idx - 1
  END DO

  IF (bot_idx == 1) THEN
    CALL finish('mo_art_prescribed_state:art_set_bot_upp_indices',      &
           &    'bottom height greater than uppermost interface level.')
  END IF
  
  DO WHILE ((upp_idx <= bot_idx)   &
       &     .AND. (zifc_sea(upp_idx) > upper_height))
    upp_idx = upp_idx + 1
  END DO

  upp_idx = upp_idx - 1


  ! compute global maximum bottom index and global minimum upper index 
  ! to ensure that they are used uniformly all over the globe
  bot_idx_max = p_max(bot_idx, p_comm_work)
  bot_idx = bot_idx_max

  upp_idx_min = p_min(upp_idx, p_comm_work)
  upp_idx = upp_idx_min

  IF (bot_idx < upp_idx) THEN
    CALL finish('mo_art_prescribed_state:art_set_bot_upp_indices',      &
           &    'bottom height smaller than upper height.')
  END IF

  
END SUBROUTINE art_set_bot_upp_indices
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_prescribe_tracers(tracer,prescr_list,p_patch,current_date)
!<
! SUBROUTINE art_prescribe_tracers
! This subroutine performs the prescription of the tracers
! Part of Module: mo_art_prescribed_state
! Author: Michael Weimer, KIT
! Initial Release: 2017-09-06
! Modifications:
!>
  IMPLICIT NONE
  REAL(wp), INTENT(INOUT) :: &
    &  tracer(:,:,:,:)        !< tracer mass mixing ratios (kg/kg)
  TYPE(t_art_prescr_list), INTENT(IN) :: &
    &  prescr_list            !< linked list with the dataset meta data
  TYPE(t_patch), INTENT(IN) :: &
    &  p_patch                !< patch on which computation is performed
  TYPE(datetime), INTENT(IN), POINTER :: &
    &  current_date           !< current simulation time
  ! local variables
  TYPE(t_art_prescr_list_element), POINTER :: &
    &  prescr_element         !< element of the list
  INTEGER :: &
    &  ijsp                   !< loop index
  INTEGER ::                &
    &  jc,jk,jb, i_startidx,& !< loop indices
    &  i_endidx,            &
    &  jg
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo

  IF (prescr_list%num_elements > 0) THEN
    jg = p_patch%id

    prescr_element => prescr_list%first_prescr_element
    art_atmo => p_art_data(jg)%atmo

    ! go through all elements of the linked list, read and vertically
    ! interpolate the data if necessary and interpolate temporally to simulation
    ! time
    DO WHILE(ASSOCIATED(prescr_element))
      IF (prescr_element%prescr%type_is_init) THEN
        CALL art_read_emissions(p_patch,current_date,prescr_element%prescr)
      ELSE
        CALL art_init_emission_struct(p_patch,current_date,                  &
                               &      prescr_element%prescr)
      END IF

      ! prescribe the tracers in the given height range

      DO ijsp = 1,prescr_element%prescr%num_vars

        DO jb = art_atmo%i_startblk, art_atmo%i_endblk
          CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

          DO jk = prescr_element%prescr%var_dep(ijsp)%upp_idx,  &
              &   prescr_element%prescr%var_dep(ijsp)%bot_idx
            DO jc = i_startidx, i_endidx
              
              tracer(jc,jk,jb,prescr_element%prescr%var_dep(ijsp)%tracer_idx) = &
                &  prescr_element%prescr%vinterp_3d(jc,jk,jb,ijsp,2)

            END DO
          END DO
        END DO

      END DO

      prescr_element => prescr_element%next_prescr_element
    END DO
  END IF
END SUBROUTINE art_prescribe_tracers
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_delete_prescr_list(this_list)
!<
! SUBROUTINE art_delete_prescr_list
! This subroutine deletes the linked list with meta data 
! of the dataset for prescribing tracers
! Part of Module: mo_art_prescribed_state
! Author: Michael Weimer, KIT
! Initial Release: 2020-06-16
! Modifications:
!>
  IMPLICIT NONE
  TYPE(t_art_prescr_list), INTENT(INOUT), TARGET :: &
    &  this_list                 !< linked list to be created
  ! local variables
  TYPE(t_art_prescr_list_element), POINTER :: &
    &  prescr_element,                        &  !< element of the list
    &  temp_element


  ! Delete all elements from 2 to end
  IF (this_list%num_elements > 1) THEN
    DO WHILE (ASSOCIATED(this_list%first_prescr_element%next_prescr_element))
      prescr_element => this_list%first_prescr_element%next_prescr_element

      DEALLOCATE(prescr_element%prescr%vinterp_3d)
      DEALLOCATE(prescr_element%prescr%var_dep)

      CALL deallocate_chem_init_chem(prescr_element%prescr%prescribed_state)
      CALL deallocate_chem_init_atm(prescr_element%prescr%prescribed_state)

      
      this_list%first_prescr_element%next_prescr_element => prescr_element%next_prescr_element
      DEALLOCATE(prescr_element)
      this_list%num_elements = this_list%num_elements - 1
    END DO
  END IF

  ! Delete the first element
  IF (this_list%num_elements == 1) THEN
    prescr_element => this_list%first_prescr_element

    DEALLOCATE(prescr_element%prescr%vinterp_3d)
    DEALLOCATE(prescr_element%prescr%var_dep)

    CALL deallocate_chem_init_chem(prescr_element%prescr%prescribed_state)
    CALL deallocate_chem_init_atm(prescr_element%prescr%prescribed_state)

    
    DEALLOCATE(prescr_element)
    NULLIFY(this_list%first_prescr_element)
    this_list%num_elements = this_list%num_elements - 1
  END IF

END SUBROUTINE art_delete_prescr_list
END MODULE mo_art_prescribed_state
