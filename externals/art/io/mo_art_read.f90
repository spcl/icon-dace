!
! mo_art_read
! This module provides read routines for 2D and 3D fields.
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

MODULE mo_art_read
! ICON
  USE mo_mpi,                           ONLY: my_process_is_mpi_workroot, p_bcast,  &
    &                                         process_mpi_root_id, p_comm_work
  USE mo_kind,                          ONLY: wp
  USE mo_model_domain,                  ONLY: p_patch
  USE mo_util_cdi,                      ONLY: read_cdi_2d, read_cdi_3d,               &
    &                                         t_inputparameters, makeinputparameters, &
    &                                         deleteInputParameters,                  &
    &                                         test_cdi_varid, get_cdi_varid
  USE mo_util_uuid_types,               ONLY: t_uuid, uuid_string_length
  USE mo_util_uuid,                     ONLY: OPERATOR(==), uuid_unparse
  USE mo_exception,                     ONLY: message, message_text, finish
! JF:   USE mo_initicon_types,                ONLY: ana_varnames_dict
  USE mo_run_config,                    ONLY: check_uuid_gracefully
  USE mo_initicon_config,               ONLY: nml_filetype => filetype
  USE mo_cdi,                           ONLY: FILETYPE_NC2, FILETYPE_NC4, FILETYPE_GRB2, &
    &                                         streaminqvlist, vlistinqvargrid,           &
    &                                         streamopenread, streamclose, gridInqUUID
  USE mo_io_util,                       ONLY: get_filetype
  USE mo_dictionary,                    ONLY: DICT_MAX_STRLEN
! ART

  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: art_read
  PUBLIC :: art_open_cdi
  PUBLIC :: art_close_cdi
  
INTERFACE art_read
  MODULE PROCEDURE art_read_2d_real
  MODULE PROCEDURE art_read_2d_int
  MODULE PROCEDURE art_read_3d_real
END INTERFACE art_read
    
CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_read_2d_real(jg, cdi_artdataset_id, cdi_param, varname, field_2d_real)
!<
! SUBROUTINE art_read_2d_real
! This subroutine reads 2d real fields from a given cdi file
! Part of Module: mo_art_read
! Author: Daniel Rieger, KIT
! Initial Release: 2014-05-12
! Modifications:
! 2015-02-11: Daniel Rieger, KIT
! - Check for availability of variable in file and print warning if not 
!   available (instead of an error)
! 2015-08-13: Daniel Rieger, KIT
! - Check for UUID of grid
! 2017-05-22: Michael Weimer, KIT
! - things related to cdi_artdataset_id have to be performed only on workroot
!   communicator. Reading the file on all processes (because data is broadcast
!   within cdi_read_3d to all processes).
!>
  INTEGER, INTENT(in)         :: &
    &  jg
  INTEGER,INTENT(in)          :: &
    &  cdi_artdataset_id           !< CDI ID of the dataset where the variable can be found
  TYPE(t_inputParameters),INTENT(inout) :: &
    &  cdi_param                   !< Parameters for read_cdi call
  CHARACTER(LEN=*),INTENT(in) :: &
    &  varname                     !< Name of the variable to be read
  REAL(wp),INTENT(inout)      :: &
    &  field_2d_real(:,:)          !< Array to store the data which is read
  !Local Variables
  INTEGER           :: &
    &  vlist_id,       & !< Variable list ID in ART dataset
    &  varid,          & !< Variable ID in ART dataset
    &  cdigridid         !< Grid ID in ART dataset
  CHARACTER(len=uuid_string_length) :: &
    &  grid_uuid_unparsed,             &  !< unparsed grid uuid (human readable)
    &  artdataset_uuid_unparsed           !< same for art-file uuid
  TYPE(t_uuid)      :: &
    &  artdataset_uuidOfHGrid             !< same, but converted to TYPE(t_uuid)
  LOGICAL           :: &
    &  lmatch,         &   !< for comparing UUIDs
    &  lvarname_in_file

  lvarname_in_file = .FALSE.
  
  IF (my_process_is_mpi_workroot()) THEN
    IF (test_cdi_varID(cdi_artdataset_id,TRIM(varname)) /= -1) THEN
      ! Check if variable exists in file
      ! Compare UUID of ART input file with UUID of grid.
      lvarname_in_file = .TRUE.
      varid     = get_cdi_varid(cdi_artdataset_id, TRIM(varname))
      vlist_id  = streaminqvlist(cdi_artdataset_id)
      cdigridid = vlistinqvargrid(vlist_id, varid)
      CALL gridInqUUID(cdigridid, artdataset_uuidOfHGrid%DATA)
      !
      ! --- compare UUID of horizontal grid file with UUID from extpar file
      lmatch = (p_patch(jg)%grid_uuid == artdataset_uuidOfHGrid)
      IF (.NOT. lmatch) THEN
        CALL uuid_unparse(p_patch(jg)%grid_uuid, grid_uuid_unparsed)
        CALL uuid_unparse(artdataset_uuidOfHGrid, artdataset_uuid_unparsed)
        WRITE(message_text,'(5a)') 'UUIDs do not match:: ICON: ', TRIM(grid_uuid_unparsed), &
          &                        ', ART: ', TRIM(artdataset_uuid_unparsed), '.'
        IF (check_uuid_gracefully) THEN
          CALL message('mo_art_read: art_read_2d_real', TRIM(message_text))
        ELSE
          CALL finish('mo_art_read: art_read_2d_real', TRIM(message_text))
        ENDIF
      ENDIF
    ELSE ! Variable is not contained in file
      CALL message('','WARNING: Variable '//TRIM(varname)//' not found in file.')
    ENDIF
  END IF
  
  ! Communicate lvarname_in_file from workroot to all processes
  ! so that only a warning is printed if variable is not found in the file
  ! (which can only be checked by the workroot communicator)
  CALL p_bcast(lvarname_in_file,process_mpi_root_id,p_comm_work)

  IF (lvarname_in_file) THEN
    CALL read_cdi_2d(cdi_param, TRIM(varname), field_2d_real)
  END IF
END SUBROUTINE art_read_2d_real
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_read_2d_int(jg, cdi_artdataset_id, cdi_param, varname, field_2d_int)
!<
! SUBROUTINE art_read_2d_int
! This subroutine reads 2d integer fields from a given cdi file
! Part of Module: mo_art_read
! Author: Daniel Rieger, KIT
! Initial Release: 2014-05-12
! Modifications:
! 2015-02-11: Daniel Rieger, KIT
! - Check for availability of variable in file and print warning if not 
!   available (instead of an error)
! 2015-08-13: Daniel Rieger, KIT
! - Check for UUID of grid
! 2017-05-22: Michael Weimer, KIT
! - things related to cdi_artdataset_id have to be performed only on workroot
!   communicator. Reading the file on all processes (because data is broadcast
!   within cdi_read_3d to all processes).
!>
  INTEGER, INTENT(in)         :: &
    &  jg
  INTEGER,INTENT(in)          :: &
    &  cdi_artdataset_id           !< CDI ID of the dataset where the variable can be found
  TYPE(t_inputparameters),INTENT(inout) :: &
    &  cdi_param                   !< Parameters for read_cdi call
  CHARACTER(LEN=*),INTENT(in) :: &
    &  varname                     !< Name of the variable to be read
  INTEGER,INTENT(inout)       :: &
    &  field_2d_int(:,:)           !< Array to store the data which is read
  !Local Variables
  INTEGER           :: &
    &  vlist_id,       & !< Variable list ID in ART dataset
    &  varid,          & !< Variable ID in ART dataset
    &  cdigridid         !< Grid ID in ART dataset
  CHARACTER(len=uuid_string_length) :: &
    &  grid_uuid_unparsed,             &  !< unparsed grid uuid (human readable)
    &  artdataset_uuid_unparsed           !< same for art-file uuid
  TYPE(t_uuid)      :: &
    &  artdataset_uuidofhgrid             !< same, but converted to TYPE(t_uuid)
  LOGICAL           :: &
    &  lmatch,         &   !< for comparing UUIDs
    &  lvarname_in_file

  lvarname_in_file = .FALSE.
  
  IF (my_process_is_mpi_workroot()) THEN
    IF (test_cdi_varid(cdi_artdataset_id,TRIM(varname)) /= -1) THEN
      ! Check if variable exists in file
      ! Compare UUID of ART input file with UUID of grid.
      lvarname_in_file = .TRUE.

      varid     = get_cdi_varid(cdi_artdataset_id, TRIM(varname))
      vlist_id  = streaminqvlist(cdi_artdataset_id)
      cdigridid = vlistinqvargrid(vlist_id, varid)
      CALL gridInqUUID(cdigridid, artdataset_uuidOfHGrid%DATA)
      
      ! --- compare UUID of horizontal grid file with UUID from extpar file
      lmatch = (p_patch(jg)%grid_uuid == artdataset_uuidofhgrid)
      IF (.NOT. lmatch) THEN
        CALL uuid_unparse(p_patch(jg)%grid_uuid, grid_uuid_unparsed)
        CALL uuid_unparse(artdataset_uuidofhgrid, artdataset_uuid_unparsed)
        WRITE(message_text,'(5a)') 'UUIDs do not match:: ICON: ', TRIM(grid_uuid_unparsed), &
          &                        ', ART: ', TRIM(artdataset_uuid_unparsed), '.'
        IF (check_uuid_gracefully) THEN
          CALL message('mo_art_read: art_read_2d_int', TRIM(message_text))
        ELSE
          CALL finish('mo_art_read: art_read_2d_int', TRIM(message_text))
        ENDIF
      ENDIF
    ELSE ! Variable is not contained in file
      CALL message('','WARNING: Variable '//TRIM(varname)//' not found in file.')
    ENDIF
  END IF

  ! Communicate lvarname_in_file from workroot to all processes
  ! so that only a warning is printed if variable is not found in the file
  ! (which can only be checked by the workroot communicator)
  CALL p_bcast(lvarname_in_file,process_mpi_root_id,p_comm_work)

  IF (lvarname_in_file) THEN
    CALL read_cdi_2d(cdi_param, TRIM(varname), field_2d_int)
  END IF
  
END SUBROUTINE art_read_2d_int
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_read_3d_real(jg, cdi_artdataset_id, cdi_param, varname, field_3d_real, &
  &                         nlev, cdi_filetype)
!<
! SUBROUTINE art_read_3d_real
! This subroutine reads 3d real fields from a given cdi file
! Part of Module: mo_art_read
! Author: Daniel Rieger, KIT
! Initial Release: 2014-05-12
! Modifications:
! 2015-02-11: Daniel Rieger, KIT
! - Check for availability of variable in file and print warning if not 
!   available (instead of an error)
! 2015-08-13: Daniel Rieger, KIT
! - Check for UUID of grid
! 2017-05-22: Michael Weimer, KIT
! - things related to cdi_artdataset_id have to be performed only on workroot
!   communicator. Reading the file on all processes (because data is broadcast
!   within cdi_read_3d to all processes).
!>
  INTEGER, INTENT(in)         :: &
    &  jg
  INTEGER,INTENT(in)          :: &
    &  cdi_artdataset_id           !< CDI ID of the dataset where the variable can be found
  TYPE(t_inputParameters),INTENT(inout) :: &
    &  cdi_param                   !< Parameters for read_cdi call
  CHARACTER(LEN=*),INTENT(in) :: &
    &  varname                     !< Name of the variable to be read
  REAL(wp),INTENT(inout)      :: &
    &  field_3d_real(:,:,:)        !< Array to store the data which is read
  INTEGER,INTENT(in)          :: &
    &  nlev                        !< number of vertical levels in dataset
  INTEGER, INTENT(in), OPTIONAL :: &
    &  cdi_filetype                !< One of CDI's FILETYPE_XXX constants
  !Local Variables
  INTEGER           :: &
    &  vlist_id,       & !< Variable list ID in ART dataset
    &  varid,          & !< Variable ID in ART dataset
    &  cdigridid         !< Grid ID in ART dataset
  CHARACTER(LEN=DICT_MAX_STRLEN)  :: mapped_name
  CHARACTER(len=uuid_string_length) :: &
    &  grid_uuid_unparsed,             &  !< unparsed grid uuid (human readable)
    &  artdataset_uuid_unparsed           !< same for art-file uuid
  TYPE(t_uuid)      :: &
    &  artdataset_uuidofhgrid             !< same, but converted to TYPE(t_uuid)
  LOGICAL           :: &
    &  lmatch,         &  !< for comparing UUIDs
    &  lvarname_in_file

  lvarname_in_file = .FALSE.

  IF (PRESENT(cdi_filetype)) THEN
    SELECT CASE(cdi_filetype)
    CASE (FILETYPE_NC2, FILETYPE_NC4)
      ! Trivial name mapping for NetCDF file
      mapped_name = TRIM(varname)
    CASE (FILETYPE_GRB2)
      ! Search name mapping for name in GRIB2 file
! JF: ! The dictionary 'ana_varnames_dict' is not yet initialized...
! JF: ! TODO: CALL art_init_aero([...]) after CALL initVarnamesDict(ana_varnames_dict)
! JF: mapped_name = TRIM(ana_varnames_dict%get(varname, default=varname))
      SELECT CASE(TRIM(varname))
      CASE('pollalnu')
        mapped_name = 'ALNUsnc'
      CASE('pollcory')
        mapped_name = 'CORYsnc'
      CASE('pollambr')
        mapped_name = 'AMBRsnc'
      CASE('pollbetu')
        mapped_name = 'BETUsnc'
      CASE('pollpoac')
        mapped_name = 'POACsnc'
      CASE('soot')
        mapped_name = 'bcarb'
      CASE('soot0')
        mapped_name = 'bcarb0'
      CASE DEFAULT
        mapped_name = TRIM(varname)
      END SELECT
    CASE DEFAULT
      CALL finish('mo_art_read: art_read_3d_real', 'Unknown file type')
    END SELECT
  ELSE
    mapped_name = TRIM(varname)
  END IF

  ! Check if variable exists in file
  ! Compare UUID of ART input file with UUID of grid.
  IF (my_process_is_mpi_workroot()) THEN
    IF (test_cdi_varid(cdi_artdataset_id,TRIM(mapped_name)) /= -1) THEN
      lvarname_in_file = .TRUE.

      varid     = get_cdi_varid(cdi_artdataset_id, TRIM(mapped_name))
      vlist_id  = streaminqvlist(cdi_artdataset_id)
      cdigridid = vlistinqvargrid(vlist_id, varid)
      CALL gridInqUUID(cdigridid, artdataset_uuidOfHGrid%DATA)
      
      ! --- compare UUID of horizontal grid file with UUID from extpar file
      lmatch = (p_patch(jg)%grid_uuid == artdataset_uuidofhgrid)
      IF (.NOT. lmatch) THEN
        CALL uuid_unparse(p_patch(jg)%grid_uuid,  grid_uuid_unparsed)
        CALL uuid_unparse(artdataset_uuidofhgrid, artdataset_uuid_unparsed)
        WRITE(message_text,'(5a)') 'UUIDs do not match:: ICON: ', TRIM(grid_uuid_unparsed), &
          &                        ', ART: ', TRIM(artdataset_uuid_unparsed), '.'
        IF (check_uuid_gracefully) THEN
          CALL message('mo_art_read: art_read_3d_real', TRIM(message_text))
        ELSE
          CALL finish('mo_art_read: art_read_3d_real', TRIM(message_text))
        ENDIF
      ENDIF
    ELSE ! Variable is not contained in file
      CALL message('','WARNING: Variable '//TRIM(mapped_name)//' not found in file.')
    ENDIF
  END IF


  ! Communicate lvarname_in_file from workroot to all processes
  ! so that only a warning is printed if variable is not found in the file
  ! (which can only be checked by the workroot communicator)
  CALL p_bcast(lvarname_in_file,process_mpi_root_id,p_comm_work)

  ! Perform CDI read operation
  IF (lvarname_in_file) THEN
    CALL read_cdi_3d(cdi_param, TRIM(mapped_name),nlev, field_3d_real)
  END IF
  
END SUBROUTINE art_read_3d_real
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_open_cdi(jg,dataset,cdi_artdataset_id,cdi_param,cdi_filetype)
!<
! SUBROUTINE art_open_cdi
! Opens a CDI dataset and creates the parameters required for reading
! this dataset. This wrapper is intended to separate ART from CDI 
! specific code and use statements.
! Part of Module: mo_art_read
! Author: Daniel Rieger, KIT
! Initial Release: 2015-02-11
! Modifications:
! 2017-05-22: Michael Weimer, KIT
! - opening file only on workroot communicator
!>
  INTEGER, INTENT(in)               :: &
    &  jg
  CHARACTER(LEN=*),INTENT(IN)       :: &
    &  dataset                           !< Path+filename of dataset for initialization of aerosol
  INTEGER, INTENT(OUT)              :: &
    &  cdi_artdataset_id                 !< CDI Dataset ID; Reference ID for input from this file
  TYPE(t_inputParameters), INTENT(OUT) :: &
    &  cdi_param                         !< Parameters for read_cdi call
  INTEGER, INTENT(OUT), OPTIONAL    :: &
    &  cdi_filetype                      !< One of CDI's FILETYPE_XXX constants
  INTEGER                           :: &
    &  index_point                       !< index of point in grid_part to exclude suffix
  CHARACTER(LEN=4)                  :: &
    &  suffix                            !< Filename suffix
  LOGICAL                           :: &
    &  l_exist                           !< Return value whether dataset exists at the chosen path
  
  CALL message('','ART: Reading data from: '//TRIM(dataset))

  cdi_artdataset_id = -1
  
  suffix=''
  
  index_point = INDEX(dataset, '.nc', back=.TRUE.)
  
  IF (index_point > 0) THEN
    ! NetCDF file
    INQUIRE (FILE=TRIM(dataset), EXIST=l_exist)
  ELSE
    ! unclear filetype
    INQUIRE (FILE=TRIM(dataset), EXIST=l_exist)
    IF (.NOT.l_exist) THEN
      INQUIRE (FILE=TRIM(dataset)//'.nc', EXIST=l_exist)
      IF (.NOT.l_exist) THEN
        INQUIRE (FILE=TRIM(dataset)//'.grb', EXIST=l_exist)
        IF (l_exist) suffix='.grb'
      ELSE
        suffix='.nc'
      END IF
    END IF
  END IF
  
  IF (.NOT.l_exist)  CALL finish('mo_art_read: art_open_cdi','Dataset '//TRIM(dataset) &
                      &          //' not found.')
  
  IF (PRESENT(cdi_filetype)) THEN
    IF (nml_filetype == -1) THEN
      cdi_filetype = get_filetype(TRIM(dataset)//TRIM(suffix)) ! determine filetype
    ELSE
      cdi_filetype = nml_filetype
    END IF
  END IF
  
  IF (my_process_is_mpi_workroot()) THEN
    ! GET THE CDI DATASET ID
    cdi_artdataset_id = streamopenread(TRIM(dataset)//TRIM(suffix))
  END IF

  CALL p_bcast(cdi_artdataset_id, process_mpi_root_id, p_comm_work)
  cdi_param = makeinputparameters(cdi_artdataset_id, p_patch(jg)%n_patch_cells_g,  &
                    &             p_patch(jg)%comm_pat_scatter_c)
  
END SUBROUTINE art_open_cdi
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_close_cdi(cdi_artdataset_id, cdi_param)
!<
! SUBROUTINE art_close_cdi
! Closes a CDI dataset. This wrapper is intended to separate ART 
! from CDI specific code and use statements.
! Part of Module: mo_art_read
! Author: Daniel Rieger, KIT
! Initial Release: 2015-02-11
! Modifications:
! 2017-05-22: Michael Weimer, KIT
! - added deletion of cdi_param. streamclose only for workroot
!>
  INTEGER, INTENT(IN)               :: &
    &  cdi_artdataset_id                 !< CDI Dataset ID; Reference ID for input from this file
  TYPE(t_inputParameters), INTENT(INOUT) :: &
    &  cdi_param                         !< Parameters for read_cdi call
  
  IF (my_process_is_mpi_workroot()) THEN
    CALL streamclose(cdi_artdataset_id)
  END IF

  CALL deleteInputParameters(cdi_param)
  
END SUBROUTINE art_close_cdi
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_read
