!
! mo_art_read_xml
! This module reads the meta data from XML files.
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

MODULE mo_art_read_xml
  ! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_grid_config,                   ONLY: nroot
  USE mo_art_config,                    ONLY: art_config, IART_PATH_LEN
  USE mo_exception,                     ONLY: message,finish,message_text
  USE mtime,                            ONLY: newDatetime
  USE mo_key_value_store,               ONLY: t_key_value_store
  ! ART
  USE mo_art_impl_constants,            ONLY: IART_XMLTAGLEN,            &
                                          &   IART_VARNAMELEN
  USE tixi,                             ONLY: tixiOpenDocument,          &
                                          &   tixiCloseDocument,         &
                                          &   tixiGetDoubleElement,      &
                                          &   tixiGetIntegerElement,     &
                                          &   tixiGetNumberOfChilds,     &
                                          &   tixiGetNamedChildrenCount, &
                                          &   tixiGetTextElement,        &
                                          &   tixiGetChildNodename,      &
                                          &   tixiGetTextAttribute,      &
                                          &   tixiGetIntegerAttribute,   &
                                          &   tixiGetDoubleAttribute,    &
                                          &   tixiGetNumberOfAttributes, &
                                          &   tixiGetAttributeName,      &
                                          &   SUCCESS
  USE mo_art_prescribed_types,          ONLY: t_art_emiss_prescribed
  USE mo_art_io_constants,              ONLY: IRES_LEN, art_get_abbr_string,  &
                                          &   IART_FILENAMELEN,               &
                                          &   STORAGE_ATTR_SEP
  USE mo_art_create_filenames,          ONLY: art_get_res_string
  USE mo_art_config,                    ONLY: art_config

  IMPLICIT NONE

  PRIVATE

  TYPE t_xml_file
    CHARACTER(LEN=IART_PATH_LEN) :: path
    INTEGER                      :: handle
  END TYPE t_xml_file

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_read_xml'
  
  PUBLIC :: t_xml_file
  PUBLIC :: art_open_xml_file
  PUBLIC :: art_get_childnumber_xml
  PUBLIC :: art_check_tracer_children_xml
  PUBLIC :: art_close_xml_file
  PUBLIC :: art_read_visparams_xml
  PUBLIC :: art_read_dia_factors_xml
  PUBLIC :: art_get_prescr_meta
  PUBLIC :: art_read_emission_dataset
  PUBLIC :: art_read_elements_xml
  PUBLIC :: art_get_text_attribute_xml


CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_open_xml_file(filename, tixi_file)
!<
! SUBROUTINE art_open_xml_file
! This subroutines opens a XML file and returns the tixi handle
! Part of Module: mo_art_read_xml
! Author: Daniel Rieger, KIT
! Initial Release: 2016-07-29
! Modifications:
! YYYY-MM-DD: <name>,<institution>
! - <description>
!>
  
  CHARACTER(LEN=*), INTENT(in)  :: &
    &  filename                      !< Path and filename of XML file
  TYPE(t_xml_file), INTENT(out) :: &
    &  tixi_file                     !< File name and handle of XML file
! Local variables
  INTEGER            :: &
    &  retval             !< Tixi return value

  tixi_file%path = TRIM(filename)

  WRITE (message_text,*) 'ART: TIXI opening XML-File ',TRIM(tixi_file%path)
  CALL message (TRIM(routine)//':art_open_xml_file', message_text)

  retval = tixiOpenDocument(TRIM(tixi_file%path), tixi_file%handle)

  IF (retval /= 0) THEN
    CALL finish(TRIM(routine)//':art_open_xml_file',  &
           &    'Could not open XML document '//TRIM(tixi_file%path)//'.')
  ENDIF
END SUBROUTINE art_open_xml_file
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_get_childnumber_xml (tixi_file, x_path, number)
!<
! SUBROUTINE art_get_childnumber_xml
! This subroutines gets the number of childs contained in a xml file
! under a certain xpath
! Part of Module: mo_art_read_xml
! Author: Sven Werchner, KIT
! Initial Release: 2016-01-14
! Modifications:
! 2016-07-29: Daniel Rieger, KIT
! - Moved opening of files to separate subroutine
!>
  TYPE(t_xml_file), INTENT(in)  :: &
    &  tixi_file                     !< File name and handle of XML file
  CHARACTER(LEN=*),INTENT(in)   :: &
    &  x_path                        !< x-path of where to determine the number of childs
  INTEGER, INTENT(out)          :: &
    &  number                        !< Number of childs in file
! Local Variables
  INTEGER            :: &
    &  retval             !< Tixi return value

  retval = tixiGetNumberOfChilds(tixi_file%handle, TRIM(x_path), number)

  IF (retval /= SUCCESS) THEN
    CALL finish(TRIM(routine)//':art_get_childnumber_xml',   &
           &    'Error in tixiGetNumberOfChilds for '//TRIM(tixi_file%path)//'.')
  ENDIF

END SUBROUTINE art_get_childnumber_xml
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_check_tracer_children_xml(tixi_file, x_path, tracer_element, number)
!<
! SUBROUTINE art_check_tracer_children_xml
! This subroutines gets the number of childs contained in a xml file
! under a certain xpath, and checks if all children are called as
! "tracer_element"
! Part of Module: mo_art_read_xml
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-09
! Modifications:
!>
  TYPE(t_xml_file), INTENT(in)  :: &
    &  tixi_file                     !< File name and handle of XML file
  CHARACTER(LEN=*), INTENT(in)  :: &
    &  x_path                        !< x-path of where to determine the number of childs
  CHARACTER(LEN=*),  INTENT(in) :: & 
    &  tracer_element                !< name how tracer are called
  INTEGER, INTENT(out)          :: &
    &  number                        !< Number of childs in file
! Local Variables
  INTEGER            :: &
    &  retval,          & !< Tixi return value
    &  number_named       !< number of named childs

  retval = tixiGetNumberOfChilds(tixi_file%handle, TRIM(x_path), number)

  IF (retval /= SUCCESS) THEN
    CALL finish(TRIM(routine)//':art_check_tracer_children_xml',   &
           &    'Error in tixiGetNumberOfChilds for '//TRIM(tixi_file%path)//'.')
  ENDIF

  retval = tixiGetNamedChildrenCount(tixi_file%handle, TRIM(x_path), tracer_element, number_named)

  IF (retval /= SUCCESS) THEN
    CALL finish(TRIM(routine)//':art_check_tracer_children_xml',   &
           &    'Error in tixiGetNamedChildrenCount for '//TRIM(tixi_file%path)//'.')
  ENDIF

  IF (number_named /= number) THEN
    CALL finish(TRIM(routine)//':art_check_tracer_children_xml',   &
           &    'There are elements that are not called '//TRIM(tracer_element)//' in ' &
           &  //TRIM(tixi_file%path)//'.')
  END IF


END SUBROUTINE art_check_tracer_children_xml
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_read_visparams_xml (jg, mode_name, ext_params, ssa_params, asy_params,     &
    &                                ext_default, ssa_default, asy_default)
!<
! SUBROUTINE art_read_visparams_xml
! This subroutine reads metadata for optical parameters of aerosols
! Part of Module: mo_art_read_xml_metadata
! Author: Sven Werchner, KIT
! Initial Release: 2015-11-13
! Modifications:
!>
  INTEGER, INTENT(IN)             :: jg
  CHARACTER(LEN=*), INTENT(IN)    :: mode_name
  REAL(wp), POINTER, INTENT(INOUT)  :: ext_params(:,:),  &
                                  &  ssa_params(:,:),  &
                                  &  asy_params(:,:),  &
                                  &  ext_default(:,:), &
                                  &  ssa_default(:,:), &
                                  &  asy_default(:,:)
  INTEGER             :: handle, retval
  INTEGER             :: jspec, npoly, i, j
  CHARACTER(LEN=100)  :: x_path_base
  CHARACTER(LEN=3)    :: i_string, j_string
  IF(art_config(jg)%iart_nonsph > 0) THEN 
    CALL message ('mo_art_read_xmldata:art_read_visparams_xml','TIXI: opening xml-file "Meng_ICON_CMD_USE.xml"')
    retval = tixiOpenDocument(TRIM(art_config(jg)%cart_input_folder)//'/Meng_ICON_CMD_USE.xml', handle)
  ELSE
    CALL message ('mo_art_read_xmldata:art_read_visparams_xml','TIXI: opening xml-file "Mie_ICON_CMD_USE.xml"')
    retval = tixiOpenDocument(TRIM(art_config(jg)%cart_input_folder)//'/Mie_ICON_CMD_USE.xml', handle)
  ENDIF
  x_path_base = "/visparams/mode[@name='"//TRIM(mode_name)//"']/"
  
  ! Checking if jspec corresponds to length of first dimension of parameter-fields
  retval = tixiGetNumberOfChilds(handle, TRIM(x_path_base), jspec)
  retval = tixiGetIntegerAttribute(handle, TRIM(x_path_base),"num_params",npoly)
  IF (jspec /= size(ext_params, 1)) THEN
    CALL finish('mo_art_read_xml_metadata:art_read_visparams_xml','1st param-field-dim doesnt match number of XML-Entries')
  ENDIF
  IF (npoly /= size(ext_params, 2)) THEN
    CALL finish('mo_art_read_xml_metadata:art_read_visparams_xml','2nd param-field-dim doesnt match mode-attribute num_params')
  ENDIF
  
  i = 1
  DO WHILE (i <= jspec)
    WRITE(i_string,"(I2)") i
    x_path_base = "/visparams/mode[@name='"//TRIM(mode_name)//"']/waveband[@number='"//TRIM(ADJUSTL(i_string))//"']/"
    j = 1
    DO WHILE (j <= npoly)
      WRITE(j_string,"(I2)") j
      retval = tixiGetDoubleAttribute(handle, TRIM(x_path_base)//"ext","p"//TRIM(ADJUSTL(j_string)),ext_params(i,j))
      retval = tixiGetDoubleAttribute(handle, TRIM(x_path_base)//"ssa","p"//TRIM(ADJUSTL(j_string)),ssa_params(i,j))
      retval = tixiGetDoubleAttribute(handle, TRIM(x_path_base)//"asy","p"//TRIM(ADJUSTL(j_string)),asy_params(i,j))
      
      retval = tixiGetDoubleAttribute(handle, TRIM(x_path_base)//"ext","default_max",ext_default(i,1))
      retval = tixiGetDoubleAttribute(handle, TRIM(x_path_base)//"ext","default_min",ext_default(i,2))
      retval = tixiGetDoubleAttribute(handle, TRIM(x_path_base)//"ssa","default_max",ssa_default(i,1))
      retval = tixiGetDoubleAttribute(handle, TRIM(x_path_base)//"ssa","default_min",ssa_default(i,2))
      retval = tixiGetDoubleAttribute(handle, TRIM(x_path_base)//"asy","default_max",asy_default(i,1))
      retval = tixiGetDoubleAttribute(handle, TRIM(x_path_base)//"asy","default_min",asy_default(i,2))     
      j = j + 1
    ENDDO
    i = i + 1
  ENDDO
  
END SUBROUTINE art_read_visparams_xml
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_read_dia_factors_xml(jg, dia_min_factor, dia_max_factor)
!<
! SUBROUTINE art_read_dia_factors_xml
! This subroutine reads metadata for optical parameters of aerosols
! Part of Module: mo_art_read_xml_metadata
! Author: Sven Werchner, KIT
! Initial Release: 2015-11-13
! Modifications:
!>
  INTEGER, INTENT(IN)     :: jg
  REAL(wp), POINTER, INTENT(INOUT)   :: dia_min_factor(:), &
                                   &  dia_max_factor(:)   
  INTEGER                 :: handle, retval
  CHARACTER(LEN=100)      :: x_path_base
      
  x_path_base = "/visparams"
  
  IF(art_config(jg)%iart_nonsph > 0) THEN 
    CALL message ('mo_art_read_xmldata:art_read_visparams_xml','TIXI: opening xml-file "Meng_ICON_CMD_USE.xml"')
    retval = tixiOpenDocument(TRIM(art_config(jg)%cart_input_folder)//'/Meng_ICON_CMD_USE.xml', handle)
  ELSE
    CALL message ('mo_art_read_xmldata:art_read_visparams_xml','TIXI: opening xml-file "Mie_ICON_CMD_USE.xml"')
    retval = tixiOpenDocument(TRIM(art_config(jg)%cart_input_folder)//'/Mie_ICON_CMD_USE.xml', handle)
  ENDIF
 
  retval = tixiGetDoubleAttribute(handle, TRIM(x_path_base),"dia_max_factor",dia_max_factor(1))
  retval = tixiGetDoubleAttribute(handle, TRIM(x_path_base),"dia_min_factor",dia_min_factor(1))
  !retval = tixiGetDoubleAttribute(handle, "/visparams", "dia_max_factor", dia_max_factor)
  !retval = tixiGetDoubleAttribute(handle, "/visparams", "dia_min_factor", dia_min_factor)
  
END SUBROUTINE art_read_dia_factors_xml
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_close_xml_file(tixi_file)
!<
! SUBROUTINE art_close_xml_file
! This subroutines closes a XML file based on the tixi handle
! Part of Module: mo_art_read_xml
! Author: Daniel Rieger, KIT
! Initial Release: 2016-07-29
! Modifications:
! YYYY-MM-DD: <name>,<institution>
! - <description>
!>
  TYPE(t_xml_file), INTENT(in)  :: &
    &  tixi_file                     !< File name and handle of XML file
! Local variables
  INTEGER            :: &
    &  retval             !< Tixi return value

  WRITE (message_text,*) 'ART: TIXI closing XML-File '//TRIM(tixi_file%path)
  CALL message (TRIM(routine)//':art_close_xml_file', message_text)

  retval = tixiCloseDocument(tixi_file%handle)
  
  IF (retval /= 0) THEN
    CALL finish(TRIM(routine)//':art_close_xml_file',      &
           &    'Could not close XML document '//TRIM(tixi_file%path)//'.')
  ENDIF

END SUBROUTINE art_close_xml_file
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_read_emission_dataset(jg,emiss_t,tixi_file_emiss)
!<
! SUBROUTINE art_read_emission_dataset
! This subroutine reads the metadata of the emission datasets given in the emnissions' xml.
! Part of Module: mo_art_read_xml
! Author: Michael Weimer, KIT
! Initial Release: 2016-10-26
! Modifications:
! yyyy-mm-dd:
! - 
!>
  IMPLICIT NONE
  INTEGER, INTENT(in)  :: &
     &  jg                                       !< p_patch%id
  TYPE(t_art_emiss_prescribed), INTENT(inout)  ::   &
     &  emiss_t                                  !< current prescribed emission dataset to be read
  TYPE(t_xml_file), INTENT(in)  ::  &
     &  tixi_file_emiss                          !< XML file with emission dataset metadata
  !local variables
  INTEGER ::              &
     &  len_io_suffix,    &
     &  len_ires_nroot,   &
     &  idxvar,           &                      !< loop index
     &  retval,           &                      !< return value of the tixi routines
     &  ival_tixi                                !< Integer return value from tixiGetIntegerElement
  CHARACTER(LEN=IART_PATH_LEN) ::  &
    &  x_path_base ,               &             !< xml internal base path (e.g. root element)
    &  x_path                                    !< xml internal path
  CHARACTER, POINTER                       :: &
    &  cval_tixi(:)                              !< Temporary storage for tracer name 
                                                 !  (required by tixi)
  CHARACTER(LEN=100)                       :: &
    &  cval,                                  &  !< Restorage of cval_tixi
    &  inventory, species                        !< temporal storage of the xml reads
  CHARACTER(LEN = 4) ::   &
    &  idxvar_str                                !< string of the current variable index 
  CHARACTER(LEN=IRES_LEN)   :: & 
    &  resolution_xml                 !< resolution attribute of the emissions root element 
                                      !  in the emission XML (!)
  CHARACTER(LEN=IART_VARNAMELEN) :: &
    &  io_suffix
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! resolution check for global domain (jg = 1) 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  x_path_base = "/emissions/"

  IF (jg == 1) THEN
    retval = tixiGetTextAttribute(tixi_file_emiss%handle,TRIM(x_path_base), &
                       &          'global_resolution',cval_tixi)

    WRITE(cval,*) cval_tixi
    IF (nroot >= 10) THEN
      len_ires_nroot = IRES_LEN+1
    ELSE
      len_ires_nroot = IRES_LEN
    END IF
    resolution_xml = cval(2:len_ires_nroot)
    len_io_suffix = len_ires_nroot+2 + LEN_TRIM(art_config(jg)%cart_io_suffix)
    io_suffix = cval(len_ires_nroot+2:len_io_suffix)

    IF ((retval /= SUCCESS) .OR. (TRIM(art_get_res_string(jg)) /= TRIM(resolution_xml))  &
    &                       .OR. (TRIM(io_suffix) /= TRIM(art_config(jg)%cart_io_suffix)) ) THEN
      CALL finish('mo_art_read_xml:art_read_emission_datasets',&
                & 'global_resolution attribute could not be read or has a wrong value: ' &
                & //TRIM(resolution_xml)//'_'//TRIM(io_suffix)//' vs. '                        &
                & //TRIM(art_get_res_string(jg))//'_'//TRIM(art_config(jg)%cart_io_suffix))
    END IF
  END IF

 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! read metadata and put them into the emission storage
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  x_path = TRIM(x_path_base)//"dataset[@id='"//TRIM(emiss_t%id)//"']"

  ! read attribute "species" from the dataset element (this could include
  ! subfolders of base path)
  retval = tixiGetTextAttribute(tixi_file_emiss%handle,TRIM(x_path),'species',cval_tixi)
  IF (retval /= SUCCESS) THEN
    CALL finish('mo_art_read_xml:art_read_emission_datasets',  &
             &   'Could not read the species attribute from dataset with id '//TRIM(emiss_t%id) &
             & //'. Does the dataset element with this id exist in'  &
             & //TRIM(tixi_file_emiss%path)//'?')
  END IF
  WRITE(cval,*) cval_tixi
  species = ADJUSTL(cval)

  ! read the "inventory" subelement
  retval = tixiGetTextElement(tixi_file_emiss%handle,TRIM(x_path)//'/inventory',cval_tixi)
  IF (retval /= SUCCESS) THEN
    CALL finish('mo_art_read_xml:art_read_emission_datasets',  &
              & 'Could not read the inventory element from dataset with id '//TRIM(emiss_t%id)//'.')
  END IF

  WRITE(cval,*) cval_tixi
  inventory = ADJUSTL(cval)

  ! read the "start_date" subelement
  retval = tixiGetTextElement(tixi_file_emiss%handle,TRIM(x_path)//'/start_date',cval_tixi)
  IF (retval /= SUCCESS) THEN
    CALL finish('mo_art_read_xml:art_read_emission_datasets',  &
             &  'Could not read the start_date element from dataset with id '  &
             & //TRIM(emiss_t%id)//'.')
  END IF

  WRITE(cval,*) cval_tixi
  emiss_t%first_date => newDatetime(TRIM(ADJUSTL(cval)))

  ! read the "end_date" subelement
  retval = tixiGetTextElement(tixi_file_emiss%handle,TRIM(x_path)//'/end_date',cval_tixi)
  IF (retval /= SUCCESS) THEN
    CALL finish('mo_art_read_xml:art_read_emission_datasets',  &
              &  'Could not read the end_date element from dataset with id '//TRIM(emiss_t%id)//'.')
  END IF


  WRITE(cval,*) cval_tixi
  emiss_t%last_date => newDatetime(TRIM(ADJUSTL(cval)))

  ! read the variable names in the netCDF files and their attribute "num_dims" 
  retval = tixiGetNumberOfChilds(tixi_file_emiss%handle, TRIM(x_path)//'/var_names',ival_tixi)
  IF ((retval /= SUCCESS) .OR. (ival_tixi <= 0)) THEN
    CALL finish('mo_art_read_xml:art_read_emission_datasets',  &
              &  'Could not read the number of childs from the var_names element.')
  END IF
  emiss_t%num_vars = ival_tixi
  ALLOCATE(emiss_t%vname(ival_tixi),emiss_t%num_dims(ival_tixi))

  x_path = TRIM(x_path)//"/var_names"
  DO idxvar = 1, emiss_t%num_vars
    WRITE(idxvar_str,'(I4)') idxvar
    idxvar_str = ADJUSTL(idxvar_str)
    retval = tixiGetTextElement(tixi_file_emiss%handle,  &
                            &   TRIM(x_path)//"/var_name["//TRIM(idxvar_str)//"]",cval_tixi)
    IF (retval /= SUCCESS) THEN
      CALL finish('mo_art_read_xml:art_read_emission_datasets',  &
              &   'Could not read the '//TRIM(idxvar_str)  &
              & //'th variable name from dataset with id '//TRIM(emiss_t%id)//'.')
    END IF
    WRITE(cval,*) cval_tixi
    emiss_t%vname(idxvar) = ADJUSTL(cval)

    retval = tixiGetIntegerAttribute(tixi_file_emiss%handle,TRIM(x_path)//"/var_name["  &
                                & //TRIM(idxvar_str)//"]",'num_dims',ival_tixi)
    IF (retval /= SUCCESS) THEN
      CALL finish('mo_art_read_xml:art_read_emission_datasets',  &
                &  'Could not read the num_dims attribute of the '//TRIM(idxvar_str) &
                &  //'th variable name from dataset with id '//TRIM(emiss_t%id)//'.')
    END IF
    IF ((ival_tixi < 2) .OR. (ival_tixi > 3)) THEN
      CALL finish('mo_art_read_xml:art_read_emission_datasets',  &
                &  'The num_dims attribute of the '//TRIM(idxvar_str)//   &
                &  'th variable name from dataset with id '//TRIM(emiss_t%id)//' has to be 2 or 3.')
    END IF
    emiss_t%num_dims(idxvar) = ival_tixi
  END DO


  ! automatically create path to the dataset
  emiss_t%path = TRIM(art_config(jg)%cart_input_folder)//'/emissions/'//TRIM(species)//'/'  &
      &  //art_get_abbr_string(emiss_t%iType_data)//'/'//TRIM(inventory)                    &
      &  //'/'//TRIM(art_get_res_string(jg))//'_'//art_config(jg)%cart_io_suffix

END SUBROUTINE art_read_emission_dataset
!
!--------------------------------------------------------------------------------------
!
SUBROUTINE art_get_prescr_meta(tixi_file,dataset_name,var_name_in_dataset, &
                    &          first_date,last_date,vert_coord,vert_unit,  &
                    &          model_name)
!<
! SUBROUTINE art_get_prescr_meta
! This subroutine reads the metadata of the external datasets for prescribing 
! tracers given in the XML file with meta data for external datsets (not
! emissions).
! Part of Module: mo_art_read_xml
! Author: Michael Weimer, KIT
! Initial Release: 2017-09-06
! Modifications:
! yyyy-mm-dd:
! - 
!>
  IMPLICIT NONE
  TYPE(t_xml_file), INTENT(in) :: &
    &  tixi_file              !< XML file with the meta data
  CHARACTER(LEN=*), INTENT(in) :: &
    &  dataset_name,              & !< name of the dataset
    &  var_name_in_dataset          !< name of the variable in the dataset
  CHARACTER(LEN=*), INTENT(out) :: &
    &  first_date,last_date,       & !< first and last date
    &  vert_coord,vert_unit,       & !< coord file and unit of the vertical axis
    &  model_name                    !< name of the external model
  ! local
  CHARACTER(LEN = IART_PATH_LEN) :: &
    &  x_path                       !< current path within the XML file
  CHARACTER, POINTER :: &
    &  cval_tixi(:)                 !< XML output
  CHARACTER(LEN = IART_XMLTAGLEN) :: & 
    &  cval                         !< XML output
  CHARACTER(LEN=*), PARAMETER :: &
    &  routine = 'mo_art_read_xml:art_get_prescr_meta'
  CHARACTER(LEN = 4) :: &
    &  num_vars_str                 !< number of variables in the dataset
  INTEGER ::   &
    &  retval,           &          !< XML return value
    &  num_vars_dataset, &          !< number of variables in the dataset
    &  ival_tixi                    !< XML output
 
  x_path = "/external_datasets/dataset[@name='"//TRIM(dataset_name)//"']"

  ! start date
  retval = tixiGetTextElement(tixi_file%handle, &
                & TRIM(x_path)//"/start_date", cval_tixi)
  IF (retval /= SUCCESS) THEN
    CALL finish(routine, 'Could not read start_date in dataset with name '  &
          &             //TRIM(dataset_name)//' in '//TRIM(tixi_file%path))
  ELSE
    WRITE(cval,*) cval_tixi
    first_date = TRIM(ADJUSTL(cval))
  END IF
  
  ! end date
  retval = tixiGetTextElement(tixi_file%handle, &
                & TRIM(x_path)//"/end_date", cval_tixi)
  IF (retval /= SUCCESS) THEN
    CALL finish(routine, 'Could not read end_date in dataset with name '   &
           &            //TRIM(dataset_name)//' in '//TRIM(tixi_file%path))
  ELSE
    WRITE(cval,*) cval_tixi
    last_date = TRIM(ADJUSTL(cval))
  END IF
  
  ! vertical coordinate file
  retval = tixiGetTextElement(tixi_file%handle, &
                & TRIM(x_path)//"/vert_coord", cval_tixi)
  IF (retval /= SUCCESS) THEN
    vert_coord = ''
  ELSE
    WRITE(cval,*) cval_tixi
    vert_coord = TRIM(ADJUSTL(cval))
  END IF
  
  ! unit of vertical axis
  retval = tixiGetTextElement(tixi_file%handle, &
                & TRIM(x_path)//"/vert_unit", cval_tixi)
  IF (retval /= SUCCESS) THEN
    vert_unit = ''
  ELSE
    WRITE(cval,*) cval_tixi
    vert_unit = TRIM(ADJUSTL(cval))
  END IF

  ! name of the external model
  retval = tixiGetTextElement(tixi_file%handle, &
                & TRIM(x_path)//"/model_name", cval_tixi)
  IF (retval /= SUCCESS) THEN
    model_name = ''
  ELSE
    WRITE(cval,*) cval_tixi
    model_name = TRIM(ADJUSTL(cval))
  END IF

  ! look for the variable called as the given parameter var_name_in_dataset
  x_path = TRIM(x_path)//"/var_names"
  retval = tixiGetNumberOfChilds(tixi_file%handle,TRIM(x_path),ival_tixi)

  IF (retval /= SUCCESS) THEN
    CALL finish(routine, 'Could not calculate number of childs of var_names '  &
             &         //'in dataset with name '//TRIM(dataset_name)//' in '   &
             &         //TRIM(tixi_file%path))
  END IF
  
  num_vars_dataset = ival_tixi
  DO WHILE (num_vars_dataset > 0)
    WRITE(num_vars_str,'(I4)') num_vars_dataset
    
    retval = tixiGetTextElement(tixi_file%handle, &
                 & TRIM(x_path)//"/var_name["//TRIM(num_vars_str)//"]",cval_tixi)
    WRITE(cval,*) cval_tixi
    cval = ADJUSTL(cval)
    IF (TRIM(var_name_in_dataset) == TRIM(cval)) THEN
      EXIT
    END IF
    
    num_vars_dataset = num_vars_dataset -1
  END DO

  IF (num_vars_dataset == 0) THEN
    CALL finish(routine, 'could not find var_name '//TRIM(var_name_in_dataset)  &
           &           //' in dataset with name '//TRIM(dataset_name)//' in '   &
           &           //TRIM(tixi_file%path))
  END IF
  
END SUBROUTINE art_get_prescr_meta
!
!--------------------------------------------------------------------------------------
!
SUBROUTINE art_read_elements_xml(tixi_file,x_path_parent,key_value_store)
!<
! SUBROUTINE art_read_elements_xml
! <description>
! Part of Module: mo_art_read_xml
! Author: Daniel Rieger, KIT
! Initial Release: 2016-12-08
! Modifications:
! 2017-07-20: Michael Weimer, KIT
! - added treatment of attributes in the storage
!>
  TYPE(t_xml_file), INTENT(in)             :: &
    &  tixi_file                                !< File name and handle of XML file
  CHARACTER(LEN=*),INTENT(in)              :: &
    &  x_path_parent                            !< X-Path to parent node (e.g. tracer)
  TYPE(t_key_value_store), INTENT(inout)   :: &
    &  key_value_store                          !< Key-Value storage to be filled with data 
                                                !  from child elements
! Local variables
  REAL(wp)                                 :: &
    &  rval_tixi                                !< Double return value from tixiGetDoubleElement
  INTEGER                                  :: &
    &  je,                                    & !< Loop counter for elements
    &  retval,                                & !< Tixi return value (i.e. error code)
    &  nelements,                             & !< number of child elements in parent node 
                                                !  (e.g. number of metadata)
    &  nattr,                                 & !< number of attributes of the current element
    &  ival_tixi,                             & !< Integer return value from tixiGetIntegerElement
    &  counter                                  !< counter for multiple equal child-tags
  CHARACTER, POINTER                       :: &
    &  cval_tixi(:)                             !< Temporary storage for character 
                                                !  (required by tixi)
  CHARACTER(LEN=IART_XMLTAGLEN)            :: &
    &  cval,                                  & !< Restorage of cval_tixi
    &  attr_name,                             & !< name of the attribute
    &  cval_elem,                             & !< Restorage of cval_tixi for metadata
    &  parent_name,                           & !< name of parent node (e.g. tracer)
    &  element_name,                          & !< Name of element (e.g. metadata)
    &  last_elem,                             & !< Name of last element (important for tracking multiple equal tags)
    &  ccounter,                              & !< character-representation of counter for same element name
    &  cid                                      !< optional additional identifier for storage
  CHARACTER(LEN=IART_FILENAMELEN)          :: &
    &  x_path_meta

  nelements = 0

  retval = tixiGetTextAttribute(tixi_file%handle,TRIM(x_path_parent),"id",cval_tixi)
  WRITE(cval,*) cval_tixi
  parent_name = TRIM(ADJUSTL(cval))

  CALL art_get_childnumber_xml(tixi_file, x_path_parent, nelements)

  IF (nelements>0) THEN
    CALL key_value_store%init(.FALSE.)
    CALL key_value_store%put("name",TRIM(parent_name))

    ! Handling multiple equal tagnames
    last_elem = ''
    counter = 1
    ccounter = ''

    DO je = 1, nelements
      ! Get the name of the element (i.e. metadata)
      retval =  tixiGetChildNodename(tixi_file%handle,x_path_parent,je,cval_tixi)
      WRITE(cval,*) cval_tixi
      element_name = TRIM(ADJUSTL(cval))
      IF (element_name(1:8) == "#comment") CYCLE
      IF(element_name == last_elem) THEN
        counter  = counter+1
      ELSE
        counter  = 1
        last_elem = element_name
      ENDIF
          
      WRITE(ccounter,*) counter
      ccounter = TRIM(ADJUSTL(ccounter))
        

      ! Determine type of metadata (real,int,char), then read and store it. Default: Character
      x_path_meta = x_path_parent//TRIM(ADJUSTL(element_name))//'['//TRIM(ADJUSTL(ccounter))//']'//'/'
      retval      = tixiGetTextAttribute(tixi_file%handle,TRIM(x_path_meta),"type",cval_tixi)
      WRITE(cval,*) cval_tixi
      IF (retval /= SUCCESS) cval = 'char'
      ! Get id of child-element (important if there are more than one)
      retval      = tixiGetTextAttribute(tixi_file%handle,TRIM(x_path_meta),"id",cval_tixi)
      WRITE(cid,*) cval_tixi
      cid = TRIM(ADJUSTL(cid))
      IF (retval /= SUCCESS) THEN
        cid = ''
      ELSE
        cid = '_'//cid
      ENDIF
      
      IF (TRIM(ADJUSTL(cval))=='real') THEN
        retval = tixiGetDoubleElement(tixi_file%handle, TRIM(x_path_meta), rval_tixi)
        CALL key_value_store%put(TRIM(ADJUSTL(element_name))//cid,rval_tixi)
      ELSEIF (TRIM(ADJUSTL(cval))=='int') THEN
        retval = tixiGetIntegerElement(tixi_file%handle, TRIM(x_path_meta), ival_tixi)
        CALL key_value_store%put(TRIM(ADJUSTL(element_name))//cid,ival_tixi)
      ELSEIF (TRIM(ADJUSTL(cval))=='char') THEN
        retval = tixiGetTextElement(tixi_file%handle, TRIM(x_path_meta), cval_tixi)
        WRITE(cval_elem,*) cval_tixi
        CALL key_value_store%put(TRIM(ADJUSTL(element_name))//cid,TRIM(ADJUSTL(cval_elem)))
      ELSE
        CALL finish(TRIM(routine)//':art_read_elements_xml', &
          &         'Meta data type '//TRIM(ADJUSTL(cval))//' of '//TRIM(parent_name)//'/' &
          &        //TRIM(element_name)//' unknown.')
      ENDIF

      ! Handle attributes of this element (of course except for the "type"
      ! attribute). Store it in the form "<attribute><STORAGE_ATTR_SEP><element_name>" in the
      ! key_value_store

      retval = tixiGetNumberOfAttributes(tixi_file%handle, &
                           & TRIM(x_path_meta),ival_tixi)

      IF (ival_tixi > 0) THEN
        nattr = ival_tixi
        DO WHILE (nattr > 0)
          retval = tixiGetAttributeName(tixi_file%handle, &
                            & TRIM(x_path_meta),nattr,cval_tixi)
          WRITE(cval,*) cval_tixi
          attr_name = TRIM(ADJUSTL(cval))
          IF (TRIM(attr_name) /= 'type') THEN
            SELECT CASE (attr_name(1:1))
              CASE('i','I')
                retval = tixiGetIntegerAttribute(tixi_file%handle,TRIM(x_path_meta),  &
                                       &         TRIM(attr_name),ival_tixi)
                IF (retval == SUCCESS) THEN
                  CALL key_value_store%put(TRIM(attr_name)//STORAGE_ATTR_SEP     &
                            &    //TRIM(ADJUSTL(element_name)),ival_tixi)
                ELSE
                  CALL finish(TRIM(routine)//':art_read_elements_xml',  &
                          &   'Could not read attribute '//TRIM(attr_name)//' of element '  &
                          & //TRIM(ADJUSTL(element_name))//' of tracer '//TRIM(parent_name))
                END IF
              CASE('r','R')
                retval = tixiGetDoubleAttribute(tixi_file%handle,TRIM(x_path_meta),  &
                                       &        TRIM(attr_name),rval_tixi)
                IF (retval == SUCCESS) THEN
                  CALL key_value_store%put(TRIM(attr_name)//STORAGE_ATTR_SEP  &
                           &     //TRIM(ADJUSTL(element_name)),rval_tixi)
                ELSE
                  CALL finish(TRIM(routine)//':art_read_elements_xml',  &
                          &   'Could not read attribute '//TRIM(attr_name)//' of element '  &
                          & //TRIM(ADJUSTL(element_name))//' of tracer '//TRIM(parent_name))
                END IF
              CASE('c','C')
                retval = tixiGetTextAttribute(tixi_file%handle,TRIM(x_path_meta),  &
                                       &         TRIM(attr_name),cval_tixi)
                WRITE(cval_elem,*) cval_tixi
                cval_elem = TRIM(ADJUSTL(cval_elem))
                IF (retval == SUCCESS) THEN
                  CALL key_value_store%put(TRIM(attr_name)//STORAGE_ATTR_SEP    &
                           &     //TRIM(ADJUSTL(element_name)),cval_elem)
                ELSE
                  CALL finish(TRIM(routine)//':art_read_elements_xml',  &
                          &   'Could not read attribute '//TRIM(attr_name)//' of element '  &
                          & //TRIM(ADJUSTL(element_name))//' of tracer '//TRIM(parent_name))
                END IF
              CASE DEFAULT
                CALL finish(TRIM(routine)//':art_read_elements_xml',    &
                  &         'Given attribute '//TRIM(attr_name)         &
                  &       //' does not begin with i,I,r,R,c or C '      &
                  &       //' which is necessary to determine its type.')
                  
            END SELECT
          END IF  ! attr_name /= type
          nattr = nattr-1
        END DO ! while nattr > 0
      END IF ! number of attributes > 0
      ! Save current child_id (used in e.g. coagulate)
      IF(cid /= '') THEN
        CALL key_value_store%put(TRIM(ADJUSTL(element_name))//ccounter//'_id',cid)
      ENDIF
    ENDDO  ! loop over elements
  ENDIF  ! number of elements > 0

END SUBROUTINE art_read_elements_xml
!
!-------------------------------------------------------------------------
!
SUBROUTINE art_get_text_attribute_xml(tixi_file,x_path, attr_name,text_attr)
!<
! SUBROUTINE art_get_text_attribute_xml
! This subroutines gets the text attribute of an element in the already opened
! XML file
! Part of Module: mo_art_read_xml
! Author: Michael Weimer, KIT
! Initial Release: 2018-07-17
! Modifications:
!>
  TYPE(t_xml_file), INTENT(in)  :: &
    &  tixi_file                     !< File name and handle of XML file
  CHARACTER(LEN=*),INTENT(in)   :: &
    &  x_path                        !< x-path of where to determine the number of childs
  CHARACTER(LEN=*),INTENT(in)   :: &
    &  attr_name                     !< name of the attribute
  CHARACTER(LEN=*), INTENT(out)          :: &
    &  text_attr                     !< attribute content
! Local Variables
  INTEGER            :: &
    &  retval             !< Tixi return value
  CHARACTER, POINTER                       :: &
    &  cval_tixi(:)                             !< Temporary storage for character
                                                !  (required by tixi)
  CHARACTER(LEN=IART_VARNAMELEN)           :: &
    &  cval_elem                                !< Restorage of cval_tixi for metadata

  retval = tixiGetTextAttribute(tixi_file%handle, TRIM(x_path), TRIM(attr_name),cval_tixi)
  WRITE(cval_elem,*) cval_tixi
  cval_elem = TRIM(ADJUSTL(cval_elem))

  IF (retval == SUCCESS) THEN
    text_attr = TRIM(cval_elem)
  ELSE
    CALL finish('mo_art_read_xml:art_get_text_attribute_xml',   &
          &     'could not read attribute '//TRIM(attr_name)    &
          &   //' from '//TRIM(x_path)//' in '//TRIM(tixi_file%path))
  END IF
  

END SUBROUTINE art_get_text_attribute_xml
!!
!!-------------------------------------------------------------------------
!!
!
!-------------------------------------------------------------------------
!
END MODULE mo_art_read_xml
