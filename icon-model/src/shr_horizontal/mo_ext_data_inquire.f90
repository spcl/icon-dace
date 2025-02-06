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

! Tools to inquire the content of external parameter input files

MODULE mo_ext_data_inquire

  USE mo_model_domain,       ONLY: t_patch
  USE mo_impl_constants,     ONLY: inwp, io3_clim, io3_ape, MAX_CHAR_LENGTH, MODIS, &
    &                              SSTICE_ANA_CLINC, GLOBCOVER2009, GLC2000
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_util_uuid_types,    ONLY: t_uuid, uuid_string_length
  USE mo_master_config,      ONLY: getModelBaseDir
  USE mo_parallel_config,    ONLY: p_test_run
  USE mo_grid_config,        ONLY: nroot, l_scm_mode
  USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config
  USE mo_radiation_config,   ONLY: albedo_type, islope_rad, irad_o3
  USE mo_extpar_config,      ONLY: itopo, extpar_filename, generate_filename, itype_lwemiss, &
    &                              t_ext_atm_attr, t_ext_o3_attr, num_lcc
  USE mo_lnd_nwp_config,     ONLY: sstice_mode
  USE mo_run_config,         ONLY: msg_level, check_uuid_gracefully, iforcing
  USE mo_io_units,           ONLY: filename_max
  USE mo_cdi,                ONLY: CDI_GLOBAL, FILETYPE_GRB2, streamOpenRead, streamInqFileType, &
    &                              streamInqVlist, cdiStringError, vlistInqVarZaxis, zaxisInqSize, &
    &                              vlistNtsteps, vlistInqVarGrid, cdiInqAttTxt, vlistInqVarIntKey, &
    &                              gridInqUUID
  USE mo_util_cdi,           ONLY: get_cdi_varID, test_cdi_varID, has_filetype_netcdf
  USE mo_util_uuid,          ONLY: OPERATOR(==), uuid_unparse
  USE mo_netcdf_errhandler,  ONLY: nf
  USE mo_netcdf
  USE mo_mpi,                ONLY: p_io, p_comm_work_test, p_comm_work, my_process_is_mpi_workroot, &
    &                              p_bcast, my_process_is_stdio


  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_ext_data_inquire'

  PUBLIC :: inquire_external_files

CONTAINS

  !-------------------------------------------------------
  !
  ! opens the external parameter file(s) and inquires
  ! file attributes
  !
  !-------------------------------------------------------
  SUBROUTINE inquire_external_files(p_patch, ext_atm_attr, ext_o3_attr)

    TYPE(t_patch),        INTENT(IN)  :: p_patch
    TYPE(t_ext_atm_attr), INTENT(OUT) :: ext_atm_attr ! extpar file attributes
    TYPE(t_ext_o3_attr),  INTENT(OUT) :: ext_o3_attr  ! o3 file attributes

    CHARACTER(len=*), PARAMETER :: &
      routine = modname//':inquire_external_files'

!--------------------------------------------------------------------------

    ! initialize file attribute object with dummy values
    !$ACC ENTER DATA CREATE(ext_atm_attr)
    !$ACC ENTER DATA CREATE(ext_o3_attr)
    CALL ext_atm_attr%init(p_patch%id)
    CALL ext_o3_attr%init(p_patch%id)

    IF (iforcing == inwp) THEN
      !------------------------------------------------!
      ! 1. Check validity of external parameter file   !
      !------------------------------------------------!
      IF ( itopo == 1 ) THEN
        CALL inquire_extpar_file(p_patch, ext_atm_attr)

        ! Debug printput
        IF ( msg_level >= 10 ) THEN
          CALL ext_atm_attr%print_values
        ENDIF
      END IF

      !------------------------------------------------!
      ! 2. Check validity of ozone file                !
      !------------------------------------------------!
      IF ( irad_o3 == io3_clim .OR. irad_o3 == io3_ape ) THEN
        CALL inquire_o3_file(p_patch, irad_o3, ext_o3_attr)

        ! Debug printput
        IF ( msg_level >= 10 ) THEN
          CALL ext_o3_attr%print_values
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE inquire_external_files


  !-------------------------------------------------------------------------
  ! Opens ExtPar file for the atmosphere and inquires file attributes.
  ! Relevant attributes and CDI IDs are stored in an object of
  ! type t_ext_atm_attr.
  !
  !-------------------------------------------------------------------------
  SUBROUTINE inquire_extpar_file(p_patch, attr)
    TYPE(t_patch),        INTENT(IN ) :: p_patch
    TYPE(t_ext_atm_attr), INTENT(OUT) :: attr    ! file attributes

    ! local variables
    CHARACTER(len=*), PARAMETER :: routine = modname//'::inquire_extpar_file'
    INTEGER :: jg
    INTEGER :: mpi_comm, vlist_id, lu_class_fraction_id, zaxis_id, var_id
    INTEGER :: horizon_id
    LOGICAL :: l_exist

    TYPE(t_uuid)                      :: extpar_uuidOfHGrid   ! uuidOfHGrid contained in the
                                                              ! extpar file
    CHARACTER(len=uuid_string_length) :: grid_uuid_unparsed   ! unparsed grid uuid (human readable)
    CHARACTER(len=uuid_string_length) :: extpar_uuid_unparsed ! same for extpar-file uuid

    LOGICAL                 :: lmatch                         ! for comparing UUIDs
    INTEGER                 :: cdiGridID

    INTEGER :: lu_var_id, localInformationNumber
    INTEGER :: ret
    CHARACTER(len=max_char_length) :: rawdata_attr

    CHARACTER(filename_max)    :: extpar_file
    INTEGER                    :: cdi_extpar_id
    INTEGER                    :: cdi_filetype
    !
    INTEGER                    :: nhori     ! number of horizon sectors
    INTEGER                    :: nclass_lu
    INTEGER                    :: nmonths_ext
    INTEGER                    :: i_lctype
    LOGICAL                    :: is_frglac_in


    !---------------------------------------------!
    ! Check validity of external parameter file   !
    !---------------------------------------------!
    IF (my_process_is_mpi_workroot()) THEN

      jg = p_patch%id

      ! generate file name
      extpar_file = TRIM(                                            &
        &             generate_filename(extpar_filename,             &
        &                               getModelBaseDir(),           &
        &                               TRIM(p_patch%grid_filename), &
        &                               nroot,                       &
        &                               p_patch%level, p_patch%id)   &
        &               )
      CALL message(routine, "extpar_file = "//TRIM(extpar_file))

      INQUIRE (FILE=extpar_file, EXIST=l_exist)
      IF (.NOT.l_exist)  CALL finish(routine,'external data file is not found.')

      ! open file
      cdi_extpar_id = streamOpenRead(TRIM(extpar_file))
      IF (cdi_extpar_id < 0) THEN
        WRITE (message_text, '(133a)') "Cannot open external parameter file ", &
             cdiStringError(cdi_extpar_id)
        CALL finish(routine, message_text)
      END IF
      cdi_filetype = streamInqFileType(cdi_extpar_id)

      vlist_id = streamInqVlist(cdi_extpar_id)

      ! get time dimension from external data file
      nmonths_ext = vlistNtsteps(vlist_id)

      IF (islope_rad(jg) >= 2) THEN
      ! get the number of horizon sectors
        horizon_id = get_cdi_varID(cdi_extpar_id, "HORIZON")
        zaxis_id   = vlistInqVarZaxis(vlist_id, horizon_id)
        nhori      = zaxisInqSize(zaxis_id)
        WRITE(message_text,'(A,I4)')  &
          & 'Number of horizon sectors in external data file = ', nhori
        CALL message(routine, message_text)
      ELSE
        nhori = -1
      ENDIF

      ! get the number of landuse classes
      lu_class_fraction_id = get_cdi_varID(cdi_extpar_id, "LU_CLASS_FRACTION")
      zaxis_id             = vlistInqVarZaxis(vlist_id, lu_class_fraction_id)
      IF (l_scm_mode) THEN
        nclass_lu = num_lcc
      ELSE
        nclass_lu = zaxisInqSize(zaxis_id)
      ENDIF

      IF ( msg_level>10 ) THEN
        WRITE(message_text,'(A,I4)')  &
          & 'Number of land_use classes in external data file = ', nclass_lu
        CALL message(routine, message_text)

        WRITE(message_text,'(A,I4)')  &
          & 'Number of months in external data file = ', nmonths_ext
        CALL message(routine, message_text)
      ENDIF

      ! make sure that num_lcc is equal to nclass_lu. If not, then the internal
      ! land-use lookup tables and the external land-use class field are inconsistent.
      IF (nclass_lu /= num_lcc) THEN
        WRITE(message_text,'(A,I3,A,I3)')  &
          & 'Number of land-use classes in external file ', nclass_lu, &
          & ' does not match ICON-internal value num_lcc ', num_lcc
        CALL finish(routine, message_text)
      ENDIF


      ! Compare UUID of external parameter file with UUID of grid.
      !
      ! get horizontal grid UUID contained in extpar file
      ! use lu_class_fraction as sample field
      cdiGridID = vlistInqVarGrid(vlist_id, lu_class_fraction_id)
      CALL gridInqUUID(cdiGridID, extpar_uuidOfHGrid%DATA)
      !
      ! --- compare UUID of horizontal grid file with UUID from extpar file
      lmatch = (p_patch%grid_uuid == extpar_uuidOfHGrid)

      IF (.NOT. lmatch) THEN
        CALL uuid_unparse(p_patch%grid_uuid, grid_uuid_unparsed)
        CALL uuid_unparse(extpar_uuidOfHGrid, extpar_uuid_unparsed)
        WRITE(message_text,'(a,a)') 'uuidOfHgrid from gridfile: ', TRIM(grid_uuid_unparsed)
        CALL message(routine, message_text)
        WRITE(message_text,'(a,a)') 'uuidOfHgrid from extpar file: ', TRIM(extpar_uuid_unparsed)
        CALL message(routine, message_text)

        WRITE(message_text,'(a)') 'Extpar file and horizontal grid file do not match!'
        IF (check_uuid_gracefully) THEN
          CALL message(routine, message_text)
        ELSE
          CALL finish(routine, message_text)
        END IF
      ENDIF


      ! Determine which data source has been used to generate the
      ! external perameters: For NetCDF format, we check the
      ! global attribute "rawdata". For GRIB2 format we check the
      ! key "localInformationNumber".
      IF (has_filetype_netcdf(cdi_filetype)) THEN
        ret      = cdiInqAttTxt(vlist_id, CDI_GLOBAL, 'rawdata', max_char_length, rawdata_attr)
        IF (INDEX(rawdata_attr,'GLC2000') /= 0) THEN
          i_lctype = GLC2000
        ELSE IF (INDEX(rawdata_attr,'GLOBCOVER2009') /= 0) THEN
          i_lctype = GLOBCOVER2009
        ELSE
          CALL finish(routine,'Unknown landcover data source')
        ENDIF
      ELSE IF (cdi_filetype == FILETYPE_GRB2) THEN
        lu_var_id              = get_cdi_varID(cdi_extpar_id, 'LU_CLASS_FRACTION')
        localInformationNumber = vlistInqVarIntKey(vlist_id, lu_var_id, "localInformationNumber")
        SELECT CASE (localInformationNumber)
        CASE (2)  ! 2 = GLC2000
          i_lctype = GLC2000
        CASE (1)  ! 1 = ESA GLOBCOVER
          i_lctype = GLOBCOVER2009
        CASE DEFAULT
          CALL finish(routine,'Unknown landcover data source')
        END SELECT
      END IF

      ! Check whether external parameter file contains MODIS albedo-data
      IF ( albedo_type == MODIS ) THEN
        IF ( (test_cdi_varID(cdi_extpar_id, 'ALB')   == -1) .OR.    &
          &  (test_cdi_varID(cdi_extpar_id, 'ALNID') == -1) .OR.    &
          &  (test_cdi_varID(cdi_extpar_id, 'ALUVD') == -1) ) THEN
          CALL finish(routine,'MODIS albedo fields missing in '//TRIM(extpar_filename))
        ENDIF
      ENDIF

      ! Check whether external parameter file contains monthly longwave surface emissivity data
      IF ( itype_lwemiss == 2 ) THEN
        IF (test_cdi_varID(cdi_extpar_id, 'EMISS')   == -1) THEN
          CALL finish(routine,'Monthly longwave surface emissvity data missing in '//TRIM(extpar_filename))
        ENDIF
      ENDIF

      ! Check whether external parameter file contains SST climatology
      IF ( sstice_mode == SSTICE_ANA_CLINC ) THEN
        IF ( test_cdi_varID(cdi_extpar_id, 'T_SEA')  == -1 ) THEN
          CALL finish(routine,'SST climatology missing in '//TRIM(extpar_filename))
        ENDIF
      ENDIF

      IF ( atm_phy_nwp_config(jg)%icpl_aero_gscp == 3  ) THEN
        ! Check whether external parameter file contains cloud droplet number climatology
        IF ( test_cdi_varID(cdi_extpar_id, 'cdnc')  == -1 ) THEN
          CALL finish(routine,'icpl_aero_gscp=3 but cloud droplet number climatology missing in '//TRIM(extpar_filename))
        ELSE
          CALL message(routine,'Found cloud droplet number in extpar file' )
        ENDIF
      ENDIF

      ! Search for glacier fraction in Extpar file
      !
      IF (has_filetype_netcdf(cdi_filetype)) THEN
        var_id = test_cdi_varID(cdi_extpar_id,'ICE')
      ELSE IF (cdi_filetype == FILETYPE_GRB2) THEN
        var_id = test_cdi_varID(cdi_extpar_id,'FR_ICE')
      ENDIF
      IF (var_id == -1) THEN
        is_frglac_in = .FALSE.
      ELSE
        is_frglac_in = .TRUE.
      ENDIF

    ENDIF ! my_process_is_mpi_workroot()

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF
    !
    ! broadcast file attributes from I-Pe to WORK Pes
    !
    ! Note: MPI broadcast of ALLOCATABLE characters
    ! is currently not supported. This is why we do use
    ! fixed length characters for the local variables.
    CALL p_bcast(cdi_filetype,  p_io, mpi_comm)
    CALL p_bcast(extpar_file,   p_io, mpi_comm)
    CALL p_bcast(cdi_extpar_id, p_io, mpi_comm)
    CALL p_bcast(nhori,         p_io, mpi_comm)
    CALL p_bcast(nclass_lu,     p_io, mpi_comm)
    CALL p_bcast(nmonths_ext,   p_io, mpi_comm)
    CALL p_bcast(is_frglac_in,  p_io, mpi_comm)
    CALL p_bcast(i_lctype,      p_io, mpi_comm)

    ! fill file attribute object
    attr%have_inquired = .TRUE.            ! document successful file inquiry
    attr%id            = p_patch%id
    attr%cdi_filetype  = cdi_filetype
    attr%extpar_file   = TRIM(extpar_file)
    attr%cdi_extpar_id = cdi_extpar_id
    attr%nhori         = nhori
    !$ACC UPDATE DEVICE(attr%nhori) ASYNC(1)
    attr%nclass_lu     = nclass_lu
    attr%nmonths_ext   = nmonths_ext
    attr%is_frglac_in  = is_frglac_in
    attr%i_lctype      = i_lctype

  END SUBROUTINE inquire_extpar_file


  !-------------------------------------------------------------------------
  ! Opens climatological O3 file and inquires file attributes.
  ! Relevant attributes are stored in an object of type t_ext_o3_attr.
  !
  !-------------------------------------------------------------------------
  SUBROUTINE inquire_o3_file(p_patch, irad_o3, attr)

    CHARACTER(len=*), PARAMETER :: &
      routine = modname//':inquire_o3_file'

    TYPE(t_patch),       INTENT(IN)  :: p_patch
    INTEGER,             INTENT(IN)  :: irad_o3
    TYPE(t_ext_o3_attr), INTENT(OUT) :: attr

    ! local
    INTEGER :: jg
    INTEGER :: mpi_comm
    INTEGER :: no_cells
    INTEGER :: ncid, dimid

    ! necessary information when reading ozone from file
    CHARACTER(filename_max)    :: ozone_file  !< file name for reading in
    CHARACTER(max_char_length) :: levelname
    CHARACTER(max_char_length) :: cellname
    CHARACTER(max_char_length) :: o3name
    CHARACTER(max_char_length) :: o3unit
    !
    LOGICAL :: l_exist
    INTEGER :: nlev_o3
    INTEGER :: nmonths

    jg = p_patch%id

    IF(irad_o3 == io3_ape) THEN
      levelname = 'level'
      cellname  = 'ncells'
      o3name    = 'O3'
      o3unit    = 'g/g'
    ELSE ! o3_clim
      levelname = 'plev'
      cellname  = 'ncells'
      o3name    = 'O3'
      o3unit    = 'g/g' !this unit ozon will have after being read out and converted from ppmv
    ENDIF

    IF_IO : IF(my_process_is_stdio()) THEN

      WRITE(ozone_file,'(a,i2.2,a)') 'o3_icon_DOM',jg,'.nc'

      ! Note resolution assignment is done per script by symbolic links

      INQUIRE (FILE=ozone_file, EXIST=l_exist)
      IF (.NOT.l_exist) THEN
        WRITE(message_text,'(a,a,a)') 'ozone file ', TRIM(ozone_file),' is not found.'
        CALL finish(routine, message_text)
      ENDIF

      !
      ! open file
      !
      CALL nf(nf90_open(TRIM(ozone_file), NF90_NOWRITE, ncid), routine)
      CALL message(routine, 'open ozone file')

      !
      ! get number of cells
      !
      CALL nf(nf90_inq_dimid (ncid, TRIM(cellname), dimid), routine)
      CALL nf(nf90_inquire_dimension(ncid, dimid, len = no_cells), routine)
      !
      WRITE(message_text,'(a,i8)') 'number of cells are', no_cells
      CALL message(routine, message_text)

      !
      ! check the number of cells and verts
      !
      IF(p_patch%n_patch_cells_g /= no_cells) THEN
        CALL finish(routine, &
          & 'Number of patch cells and cells in ozone file do not match.')
      ENDIF

      !
      ! check the time structure
      !
      CALL nf(nf90_inq_dimid (ncid, 'time', dimid), routine)
      CALL nf(nf90_inquire_dimension(ncid, dimid, len = nmonths), routine)
      WRITE(message_text,'(A,I4)')  &
        & 'Number of months in ozone file = ', nmonths
      CALL message(routine, message_text)

      !
      ! check the vertical structure
      !
      CALL nf(nf90_inq_dimid (ncid,TRIM(levelname), dimid), routine)
      CALL nf(nf90_inquire_dimension(ncid, dimid, len = nlev_o3), routine)

      WRITE(message_text,'(A,I4)')  &
        & 'Number of pressure levels in ozone file = ', nlev_o3
      CALL message(routine, message_text)

      !
      ! close file
      !
      CALL nf(nf90_close(ncid), routine)

    END IF IF_IO

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    !
    ! broadcast file attributes from I-Pe to WORK Pes
    !
    ! Note: MPI broadcast of ALLOCATABLE characters
    ! is currently not supported. This is why we do use
    ! fixed length characters for the local variables.
    CALL p_bcast(nlev_o3,   p_io, mpi_comm)
    CALL p_bcast(nmonths,   p_io, mpi_comm)
    CALL p_bcast(levelname, p_io, mpi_comm)
    CALL p_bcast(cellname,  p_io, mpi_comm)
    CALL p_bcast(o3name,    p_io, mpi_comm)
    CALL p_bcast(o3unit,    p_io, mpi_comm)


    ! fill file attribute object
    attr%have_inquired = .TRUE.     ! document successful file inquiry
    attr%id            = jg
    attr%nlev_o3       = nlev_o3
    attr%nmonths       = nmonths
    attr%levelname     = TRIM(levelname)
    attr%cellname      = TRIM(cellname)
    attr%o3name        = TRIM(o3name)
    attr%o3unit        = TRIM(o3unit)

  END SUBROUTINE inquire_o3_file

END MODULE mo_ext_data_inquire
